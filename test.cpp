// #include "mpi.h" 
#include <upcxx/upcxx.hpp>
#include <stdio.h> 

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#define LOG(message) std::cout << myId << ": " << message << std::endl;
#define LOGV(name, value) std::cout << myId << ": " << name << " = " << value << std::endl;

int getObjectCount(std::string filename);
void readData(std::string filename, upcxx::global_ptr<double> dataVector);
void splitData(int myId, int numProcs, int totalDataSize, int* ownObjectStarts, int* ownObjectEnds);
void saveData(FILE* resultFile, double* xPositionVector, double* yPositionVector, double* zPositionVector, int totalObjectCount);

void logArray(int myId, std::string name, double* array, int iStart, int iEnd)
{
	std::cout << myId << ": " << name << ": ";
	for(int i = iStart; i < iEnd; i++)
	{
		std::cout << array[i] << ", ";
	}
	std::cout << std::endl;
}

int main(int argc, char *argv[])
{
	int myId = 0, numProcs = 2;
	std::string filename = "resultdata.txt";
	FILE* resultFile;
	int totalObjectCount = 1000, ownObjectCount, totalDataSize, ownDataSize;
    const int propertyCount = 7;
	double dt = 60;			//[s]
	double Tmax = 121;//2.6e6;	//Miesiac
	double G = 6.674e-11;

	upcxx::init();
	myId = upcxx::rank_me();
	numProcs = upcxx::rank_n();

	upcxx::global_ptr<int> totalObjectCountPtr = nullptr;

	int* ownObjectStarts = new int[numProcs];
	int* ownObjectEnds = new int[numProcs];
	int ownObjectStart, ownObjectEnd;

	if (myId == 0)
	{
		totalObjectCount = getObjectCount(filename);
		totalObjectCountPtr = upcxx::new_<int>();
		upcxx::rput(totalObjectCount, totalObjectCountPtr);
	}
	totalObjectCountPtr = upcxx::broadcast(totalObjectCountPtr, 0).wait();
	totalObjectCount = upcxx::rget(totalObjectCountPtr).wait();
    totalDataSize = totalObjectCount * propertyCount;

    upcxx::global_ptr<double> dataVector = nullptr;

	if (myId == 0)
	{
        dataVector = upcxx::new_array<double>(totalDataSize);
		readData(filename, dataVector);
	}
    dataVector = upcxx::broadcast(dataVector, 0).wait();

	//Przeslanie danych poczatkowych oraz ich podzial przez indeksy.
	splitData(myId, numProcs, totalObjectCount, ownObjectStarts, ownObjectEnds);
	ownObjectStart = ownObjectStarts[myId];
	ownObjectEnd = ownObjectEnds[myId];
	ownObjectCount = ownObjectEnd - ownObjectStart + 1;
    ownDataSize = ownObjectCount * propertyCount;

    double* xPositionVector = new double[totalObjectCount];   //m
    double* yPositionVector = new double[totalObjectCount];	//m
	double* zPositionVector = new double[totalObjectCount];	//m
    double* xVelocityVector = new double[totalObjectCount];   //m/s
    double* yVelocityVector = new double[totalObjectCount];	//m/s
	double* zVelocityVector = new double[totalObjectCount];	//m/s
    double* xAccelerationVector = new double[totalObjectCount];  //m/s2
	double* yAccelerationVector = new double[totalObjectCount];	//m/s2
	double* zAccelerationVector = new double[totalObjectCount];	//m/s2
    double* massVector = new double[totalObjectCount];
    for(int i = 0; i < totalObjectCount; i++)
    {
		xPositionVector[i] = upcxx::rget(dataVector + i*propertyCount).wait();
		yPositionVector[i] = upcxx::rget(dataVector + i*propertyCount + 1).wait();
		zPositionVector[i] = upcxx::rget(dataVector + i*propertyCount + 2).wait();
		xVelocityVector[i] = upcxx::rget(dataVector + i*propertyCount + 3).wait();
		yVelocityVector[i] = upcxx::rget(dataVector + i*propertyCount + 4).wait();
		zVelocityVector[i] = upcxx::rget(dataVector + i*propertyCount + 5).wait();
		massVector[i] = upcxx::rget(dataVector + i*propertyCount + 6).wait();
    }

	double xPosDiff;
	double yPosDiff;
	double zPosDiff;
	double r2;
	double magnitude;
	double angleH;
	double angleV;

	if(myId == 0)
	{
		resultFile = fopen(filename.c_str(), "w");
		fprintf(resultFile, "id;x;y;z\n");
	}

    upcxx::barrier();

	int writeCounter = 0;
	//Petla symulacji
	for (double t = 0; t < Tmax; t += dt, writeCounter++)
	{
		//Iteracja po cialach wlasnych danego procesu
		for (int i = ownObjectStart; i < ownObjectEnd + 1; i++)
		{
			xAccelerationVector[i] = 0;
			yAccelerationVector[i] = 0;
			zAccelerationVector[i] = 0;

			//Iteracja po wszystkich cialach symulacji
			for (int j = 0; j < totalObjectCount; j++)
			{
				if (i == j)
					continue;

				//Obliczenie przyspieszenia jakie cialo j wywiera na cialo wlasne i
				xPosDiff = xPositionVector[i] - xPositionVector[j];
				yPosDiff = yPositionVector[i] - yPositionVector[j];
				zPosDiff = zPositionVector[i] - zPositionVector[j];
				r2 = pow(xPosDiff, 2) + pow(yPosDiff, 2) + pow(zPosDiff, 2);
				magnitude = G*massVector[j] / r2;
				angleH = atan2(yPosDiff, sqrt(pow(zPosDiff, 2) + pow(xPosDiff, 2)));
				angleV = atan2(zPosDiff, xPosDiff);
				xAccelerationVector[i] += -magnitude*cos(angleH)*cos(angleV);
				yAccelerationVector[i] += -magnitude*sin(angleH);
				zAccelerationVector[i] += -magnitude*cos(angleH)*sin(angleV);
			}
		}

		// LOG("Applying acceleration. Before:");
		// logArray(myId, "xPos", xPositionVector, 0, totalObjectCount);
		// logArray(myId, "yPos", yPositionVector, 0, totalObjectCount);
		// logArray(myId, "zPos", zPositionVector, 0, totalObjectCount);
		// logArray(myId, "xVel", xVelocityVector, 0, totalObjectCount);
		// logArray(myId, "yVel", yVelocityVector, 0, totalObjectCount);
		// logArray(myId, "zVel", zVelocityVector, 0, totalObjectCount);
		//Zastosowanie obliczonych zmian predkosci i polozenia
		for (int i = ownObjectStart; i < ownObjectEnd + 1; i++)
		{
			xVelocityVector[i] += xAccelerationVector[i] * dt;
			yVelocityVector[i] += yAccelerationVector[i] * dt;
			zVelocityVector[i] += zAccelerationVector[i] * dt;
			xPositionVector[i] += xVelocityVector[i] * dt;
			yPositionVector[i] += yVelocityVector[i] * dt;
			zPositionVector[i] += zVelocityVector[i] * dt;
			upcxx::rput(xPositionVector[i], dataVector + i*propertyCount);
			upcxx::rput(yPositionVector[i], dataVector + i*propertyCount + 1);
			upcxx::rput(zPositionVector[i], dataVector + i*propertyCount + 2);
		}
		// LOG("Applying acceleration. After:");
		// logArray(myId, "xPos", xPositionVector, 0, totalObjectCount);
		// logArray(myId, "yPos", yPositionVector, 0, totalObjectCount);
		// logArray(myId, "zPos", zPositionVector, 0, totalObjectCount);
		// logArray(myId, "xVel", xVelocityVector, 0, totalObjectCount);
		// logArray(myId, "yVel", yVelocityVector, 0, totalObjectCount);
		// logArray(myId, "zVel", zVelocityVector, 0, totalObjectCount);
		
		upcxx::barrier();

		for(int i = 0; i < totalObjectCount; i++)
		{
			if(i < ownObjectStart || i > ownObjectEnd)
			{
				xPositionVector[i] = upcxx::rget(dataVector + i*propertyCount).wait();
				yPositionVector[i] = upcxx::rget(dataVector + i*propertyCount + 1).wait();
				zPositionVector[i] = upcxx::rget(dataVector + i*propertyCount + 2).wait();
			}
		}

		//Co 10 minut zapisanie danych do pliku
		if (myId == 0)// && writeCounter % 10 == 0)
			saveData(resultFile, xPositionVector, yPositionVector, zPositionVector, totalObjectCount);

		upcxx::barrier();
	}

	if(myId == 0)
		fclose(resultFile);

	delete[] ownObjectStarts;
	delete[] ownObjectEnds;

	delete[] xPositionVector;
	delete[] yPositionVector;
	delete[] zPositionVector;
	delete[] xVelocityVector;
	delete[] yVelocityVector;
	delete[] zVelocityVector;
	delete[] massVector;
	delete[] xAccelerationVector;
	delete[] yAccelerationVector;
	delete[] zAccelerationVector;

	if(myId == 0)
	{
		upcxx::delete_array(dataVector);
		upcxx::delete_(totalObjectCountPtr);
	}
	upcxx::finalize();

	return 0;
}

int getObjectCount(std::string filename)
{
	return 2;
	//Zliczenie ilosci wpisow w pliku danych
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

void splitData(int myId, int numProcs, int totalObjectCount, int* ownObjectStarts, int* ownObjectEnds)
{
	//Podzial danych miedzy procesy poprzez nadanie im pewnych zakresow indeksow.
	int baseCount = totalObjectCount / numProcs;
	int leftover = totalObjectCount % numProcs;

	ownObjectStarts[0] = 0;
	for (int i = 1; i < numProcs; i++)
	{
		ownObjectStarts[i] = ownObjectStarts[i - 1] + baseCount;
		if (leftover > 0)
		{
			ownObjectStarts[i] += 1;
			leftover--;
		}
		ownObjectEnds[i - 1] = ownObjectStarts[i] - 1;
	}
	ownObjectEnds[numProcs - 1] = totalObjectCount - 1;
}

void readData(std::string filename, upcxx::global_ptr<double> dataVector)
{
    // upcxx::rput(0.0, dataVector);
    // upcxx::rput(0.0, dataVector+1);
    // upcxx::rput(0.0, dataVector+2);
    // upcxx::rput(0.0, dataVector+3);
    // upcxx::rput(0.0, dataVector+4);
    // upcxx::rput(0.0, dataVector+5);
    // upcxx::rput(5.972e24, dataVector+6);

    // upcxx::rput(3.844e8, dataVector+7);
    // upcxx::rput(0.0, dataVector+8);
    // upcxx::rput(0.0, dataVector+9);
    // upcxx::rput(1.022e3, dataVector+10);
    // upcxx::rput(1.0, dataVector+11);
    // upcxx::rput(1.0, dataVector+12);
    // upcxx::rput(7.347e22, dataVector+13);

    //Wczytanie danych cial niebieskich z pliku
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
						dataVector[count*7] = atof(token.c_str());
						break;
					case 2:
						dataVector[count*7 + 1] = atof(token.c_str());
						break;
					case 3:
						dataVector[count*7 + 2] = atof(token.c_str());
						break;
					case 4:
						dataVector[count*7 + 3] = atof(token.c_str());
						break;
					case 5:
						dataVector[count*7 + 4] = atof(token.c_str());
						break;
					case 6:
						dataVector[count*7 + 5] = atof(token.c_str());
						break;
					}
					line.erase(0, pos + delimiter.length());
					datapos++;
				}
				dataVector[count*7 + 6] = atof(line.c_str());
			}
			++count;
		}
		datafile.close();
	}
}

void saveData(FILE* resultFile, double* xPositionVector, double* yPositionVector, double* zPositionVector, int totalObjectCount)
{
	for (int i = 0; i < totalObjectCount; i++)
	{
		fprintf(resultFile, "%d;%f;%f;%f\n", i, xPositionVector[i], yPositionVector[i], zPositionVector[i]);
	}
}