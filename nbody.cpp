#include "mpi.h" 
#include <stdio.h> 

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

int getDataSize(std::string filename);
void readData(std::string filename, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector);
void splitData(int myId, int numProcs, int totalDataSize, int* partialDataStarts, int* partialDataEnds);
void broadcastInitialData(int totalDataSize, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector);
void broadcastData(int myId, int numProcs, double* xPosVector, double* yPosVector, double* zPosVector, int* partialDataStarts, int* partialDataEnds);
void saveData(FILE* resultFile, double* xPosVector, double* yPosVector, double* zPosVector, int totalDataSize);

int main(int argc, char *argv[])
{
	int myId = 0, numProcs = 2;
	std::string filename = "resultdata.txt";
	FILE* resultFile;
	int totalDataSize = 1000, ownDataSize;
	double dt = 60;			//[s]
	double Tmax = 2.6e6;	//Miesiac
	double G = 6.674e-11;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

	int* partialDataStarts = new int[numProcs];
	int* partialDataEnds = new int[numProcs];
	int ownDataStart, ownDataEnd;

	if (myId == 0)
	{
		totalDataSize = getDataSize(filename);
	}

	MPI_Bcast(&totalDataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double* xPosVector = new double[totalDataSize];  //[m]
	double* yPosVector = new double[totalDataSize];  //[m]
	double* zPosVector = new double[totalDataSize];  //[m]
	double* xVelVector = new double[totalDataSize];  //[m/s]
	double* yVelVector = new double[totalDataSize];  //[m/s]
	double* zVelVector = new double[totalDataSize];  //[m/s]
	double* massVector = new double[totalDataSize];  //[kg]

	if (myId == 0)
	{
		readData(filename, xPosVector, yPosVector, zPosVector, xVelVector, yVelVector, zVelVector, massVector);
	}

	//Przes�anie danych pocz�tkowych oraz ich podzia� przez indeksy.
	broadcastInitialData(totalDataSize, xPosVector, yPosVector, zPosVector, xVelVector, yVelVector, zVelVector, massVector);
	splitData(myId, numProcs, totalDataSize, partialDataStarts, partialDataEnds);
	ownDataStart = partialDataStarts[myId];
	ownDataEnd = partialDataEnds[myId];
	ownDataSize = ownDataEnd - ownDataStart + 1;

	//Wektory przyspieszenia w�asnych dla danego procesu cia�. Nie przechowuj� cia� przynale�nych do innych proces�w.
	double* xAccelerationVector = new double[ownDataSize];  //m/s2
	double* yAccelerationVector = new double[ownDataSize];	//m/s2
	double* zAccelerationVector = new double[ownDataSize];	//m/s2

	double xPosDiff;
	double yPosDiff;
	double zPosDiff;
	double r2;
	double magnitude;
	double angleH;
	double angleV;

	resultFile = fopen(filename.c_str(), "w");
	fprintf(resultFile, "id;x;y;z\n");

	int writeCounter = 0;
	//P�tla symulacji
	for (double t = 0; t < Tmax; t += dt, writeCounter++)
	{
		//Iteracja po cia�ach w�asnych danego procesu
		for (int i = ownDataStart; i < ownDataEnd + 1; i++)
		{
			xAccelerationVector[i - ownDataStart] = 0;
			yAccelerationVector[i - ownDataStart] = 0;
			zAccelerationVector[i - ownDataStart] = 0;

			//Iteracja po wszystkich cia�ach symulacji
			for (int j = 0; j < totalDataSize; j++)
			{
				if (i == j)
					continue;

				//Obliczenie przyspieszenia jakie cia�o j wywiera na cia�o w�asne i
				xPosDiff = xPosVector[i] - xPosVector[j];
				yPosDiff = yPosVector[i] - yPosVector[j];
				zPosDiff = zPosVector[i] - zPosVector[j];
				r2 = pow(xPosDiff, 2) + pow(yPosDiff, 2) + pow(zPosDiff, 2);
				magnitude = G*massVector[i]*massVector[j] / r2;
				angleH = atan2(yPosDiff, sqrt(pow(zPosDiff, 2) + pow(xPosDiff, 2)));
				angleV = atan2(zPosDiff, xPosDiff);
				xAccelerationVector[i - ownDataStart] += -magnitude*cos(angleH)*cos(angleV) / massVector[i];
				yAccelerationVector[i - ownDataStart] += -magnitude*sin(angleH) / massVector[i];
				zAccelerationVector[i - ownDataStart] += -magnitude*cos(angleH)*sin(angleV) / massVector[i];
			}
		}

		//Zastosowanie obliczonych zmian pr�dko�ci i po�o�enia
		for (int i = ownDataStart; i < ownDataEnd + 1; i++)
		{
			xVelVector[i] += xAccelerationVector[i - ownDataStart] * dt;
			yVelVector[i] += yAccelerationVector[i - ownDataStart] * dt;
			zVelVector[i] += zAccelerationVector[i - ownDataStart] * dt;
			xPosVector[i] += xVelVector[i] * dt;
			yPosVector[i] += yVelVector[i] * dt;
			zPosVector[i] += zVelVector[i] * dt;
		}

		broadcastData(myId, numProcs, xPosVector, yPosVector, zPosVector, partialDataStarts, partialDataEnds);
		
		//Co 10 minut zapisanie danych do pliku
		if (myId == 0 && writeCounter % 10 == 0) {
			saveData(resultFile, xPosVector, yPosVector, zPosVector, totalDataSize);
		}
		writeCounter++;
	}

	fclose(resultFile);

	delete[] xPosVector;
	delete[] yPosVector;
	delete[] zPosVector;
	delete[] xVelVector;
	delete[] yVelVector;
	delete[] zVelVector;
	delete[] massVector;
	delete[] partialDataStarts;
	delete[] partialDataEnds;

	delete[] xAccelerationVector;
	delete[] yAccelerationVector;
	delete[] zAccelerationVector;

	MPI_Finalize();

	return 0;
}

int getDataSize(std::string filename)
{
	//Zliczenie ilo�ci wpis�w w pliku danych
	int count = 0;
	std::string line;
	std::ifstream datafile("nbodydata.txt");
	if (datafile.is_open()) {
		while (getline(datafile, line)) {
			++count;
		}
		datafile.close();
	}
	return (count - 1);
}

void splitData(int myId, int numProcs, int totalDataSize, int* partialDataStarts, int* partialDataEnds)
{
	//Podzia� danych mi�dzy procesy poprzez nadanie im pewnych zakres�w indeks�w.
	int baseCount = totalDataSize / numProcs;
	int leftover = totalDataSize%numProcs;

	partialDataStarts[0] = 0;
	for (int i = 1; i < numProcs; i++)
	{
		partialDataStarts[i] = partialDataStarts[i - 1] + baseCount;
		if (leftover > 0)
		{
			partialDataStarts[i] += 1;
			leftover--;
		}
		partialDataEnds[i - 1] = partialDataStarts[i] - 1;
	}
	partialDataEnds[numProcs - 1] = totalDataSize - 1;
}

void readData(std::string filename, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector)
{
	//Wczytanie danych cia� niebieskich z pliku
	std::string line;
	std::string delimiter = ";";
	int count = -1;
	int pos;
	std::ifstream datafile("nbodydata.txt");
	if (datafile.is_open()) {
		while (getline(datafile, line)) {
			int datapos = 0;
			if (count != -1) {
				while ((pos = line.find(delimiter)) != std::string::npos) {
					std::string token = line.substr(0, pos);
					switch (datapos) {
					case 1:
						xPosVector[count] = atof(token.c_str());
						break;
					case 2:
						yPosVector[count] = atof(token.c_str());
						break;
					case 3:
						zPosVector[count] = atof(token.c_str());
						break;
					case 4:
						xVelVector[count] = atof(token.c_str());
						break;
					case 5:
						yVelVector[count] = atof(token.c_str());
						break;
					case 6:
						zVelVector[count] = atof(token.c_str());
						break;
					}
					line.erase(0, pos + delimiter.length());
					datapos++;
				}
				massVector[count] = atof(line.c_str());
			}
			++count;
		}
		datafile.close();
	}
}

void broadcastInitialData(int totalDataSize, double* xPosVector, double* yPosVector, double* zPosVector, double* xVelVector, double* yVelVector, double* zVelVector, double* massVector)
{
	MPI_Bcast(xPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(yPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(zPosVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(xVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(yVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(zVelVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(massVector, totalDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void broadcastData(int myId, int numProcs, double* xPosVector, double* yPosVector, double* zPosVector, int* partialDataStarts, int* partialDataEnds)
{
	int partialDataSize;
	for (int i = 0; i < numProcs; i++)
	{
		partialDataSize = partialDataEnds[i] - partialDataStarts[i] + 1;
		MPI_Bcast(xPosVector + partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
		MPI_Bcast(yPosVector + partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
		MPI_Bcast(zPosVector + partialDataStarts[i], partialDataSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
	}
}

void saveData(FILE* resultFile, double* xPosVector, double* yPosVector, double* zPosVector, int totalDataSize)
{
	for (int i = 0; i < totalDataSize; i++)
	{
		fprintf(resultFile, "%d;%f;%f;%f\n", i, xPosVector[i], yPosVector[i], zPosVector[i]);
	}
}