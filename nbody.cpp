// #include "mpi.h" 
#include <upcxx/upcxx.hpp>
#include <stdio.h> 

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#define PROPERTIES 7

int getObjectCount(std::string filename);
void readData(std::string filename, upcxx::global_ptr<double> dataVector);
void splitData(int myId, int numProcs, int totalDataSize, int* ownObjectStarts, int* ownObjectEnds);
void initializeVectors(int totalObjectCount, upcxx::global_ptr<double> dataVector, double* xPositionVector, double* yPositionVector, double* zPositionVector, double* xVelocityVector, double* yVelocityVector, double* zVelocityVector, double* massVector);
void updateDataVector(int i, upcxx::global_ptr<double> dataVector, double* xPositionVector, double* yPositionVector, double* zPositionVector);
void updatePositionVectors(int ownObjectStart, int ownObjectEnd, int totalObjectCount, upcxx::global_ptr<double> dataVector, double* xPositionVector, double* yPositionVector, double* zPositionVector);
void saveData(FILE* resultFile, double* xPositionVector, double* yPositionVector, double* zPositionVector, int totalObjectCount);

int main(int argc, char *argv[])
{
    int myId = 0, numProcs = 2;
    std::string filename = "resultdata.txt";
    FILE* resultFile;
    int totalObjectCount = 1000, ownObjectCount, totalDataSize;
    double dt = 60;         //s
    double Tmax = 2.6e6;    //Miesiac
    double G = 6.674e-11;

    //Inicjalizacja UPC++, pobranie Id wlasanego i ilosci procesow
    upcxx::init();
    myId = upcxx::rank_me();
    numProcs = upcxx::rank_n();

    //Utworzenie globalnego wskaznika na liczbe obiektow symulacji
    upcxx::global_ptr<int> totalObjectCountPtr = nullptr;

    int* ownObjectStarts = new int[numProcs];
    int* ownObjectEnds = new int[numProcs];
    int ownObjectStart, ownObjectEnd;

    //Alokacja pamieci pod globalny wskaznik ilosci obiektow symulacji, nadanie mu wartosci
    if (myId == 0)
    {
        totalObjectCount = getObjectCount(filename);
        totalObjectCountPtr = upcxx::new_<int>();
        upcxx::rput(totalObjectCount, totalObjectCountPtr);
    }
    //Rozeslanie adresu wskaznika oraz pobranie wartosci z pod niego
    totalObjectCountPtr = upcxx::broadcast(totalObjectCountPtr, 0).wait();
    totalObjectCount = upcxx::rget(totalObjectCountPtr).wait();
    totalDataSize = totalObjectCount * PROPERTIES;

    //Utworzenie globalnego wskaznika na vektor danych
    upcxx::global_ptr<double> dataVector = nullptr;

    //Alokacja pamieci pod globalny wskaznik na wektor danych, wczytanie danych
    if (myId == 0)
    {
        dataVector = upcxx::new_array<double>(totalDataSize);
        readData(filename, dataVector);
    }
    //Rozeslanie adresu globalnego wskaznika na wektor danych
    dataVector = upcxx::broadcast(dataVector, 0).wait();

    //Podzial obiektow symulacji miedzy procesy na podstawie Id procesow
    splitData(myId, numProcs, totalObjectCount, ownObjectStarts, ownObjectEnds);
    ownObjectStart = ownObjectStarts[myId];
    ownObjectEnd = ownObjectEnds[myId];

    //Utworzenie wektorow danych symulacji
    double* xPositionVector = new double[totalObjectCount];    //m
    double* yPositionVector = new double[totalObjectCount];    //m
    double* zPositionVector = new double[totalObjectCount];    //m
    double* xVelocityVector = new double[totalObjectCount];    //m/s
    double* yVelocityVector = new double[totalObjectCount];    //m/s
    double* zVelocityVector = new double[totalObjectCount];    //m/s
    double* xAccelerationVector = new double[totalObjectCount];//m/s2
    double* yAccelerationVector = new double[totalObjectCount];//m/s2
    double* zAccelerationVector = new double[totalObjectCount];//m/s2
    double* massVector = new double[totalObjectCount];         //kg
    
    //Inicjalizacja wektorow danych z wektoru globalnego
    initializeVectors(totalObjectCount, dataVector, xPositionVector, yPositionVector, zPositionVector, xVelocityVector, yVelocityVector, zVelocityVector, massVector);

    //Przygotowanie pliku wyjsciowego do zapisu
    if(myId == 0)
    {
        resultFile = fopen(filename.c_str(), "w");
        fprintf(resultFile, "id;x;y;z\n");
    }

    //Synchronizacja procesow
    upcxx::barrier();

    int writeCounter = 0;
    double xPosDiff;
    double yPosDiff;
    double zPosDiff;
    double r2;
    double magnitude;
    double angleH;
    double angleV;
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

        //Zastosowanie obliczonych zmian predkosci i polozenia
        for (int i = ownObjectStart; i < ownObjectEnd + 1; i++)
        {
            xVelocityVector[i] += xAccelerationVector[i] * dt;
            yVelocityVector[i] += yAccelerationVector[i] * dt;
            zVelocityVector[i] += zAccelerationVector[i] * dt;
            xPositionVector[i] += xVelocityVector[i] * dt;
            yPositionVector[i] += yVelocityVector[i] * dt;
            zPositionVector[i] += zVelocityVector[i] * dt;

            updateDataVector(i, dataVector, xPositionVector, yPositionVector, zPositionVector);
        }

        upcxx::barrier();

        //Aktualizacja wektorow polozenia wszyzstkich procesow na podstawie dataVector
        updatePositionVectors(ownObjectStart, ownObjectEnd, totalObjectCount, dataVector, xPositionVector, yPositionVector, zPositionVector);

        upcxx::barrier();

        //Co 10 minut zapisanie danych do pliku
        if (myId == 0 && writeCounter % 10 == 0)
            saveData(resultFile, xPositionVector, yPositionVector, zPositionVector, totalObjectCount);
        //Co 100 minut informacja o postepie symulacji na konsole
        if (myId == 0 && writeCounter % 100 == 0)
            std::cout << t << std::endl;
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

    //Zwolnienie pamici globalnej oraz zakonczenie dzialania z UPC++
    if(myId == 0)
    {
        upcxx::delete_array(dataVector);
        upcxx::delete_(totalObjectCountPtr);
    }
    upcxx::finalize();

    return 0;
}

//Zliczenie ilosci wpisow w pliku danych
int getObjectCount(std::string filename)
{
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

//Podzial danych miedzy procesy poprzez nadanie im pewnych zakresow indeksow
void splitData(int myId, int numProcs, int totalObjectCount, int* ownObjectStarts, int* ownObjectEnds)
{
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

//Wczytanie danych cial niebieskich z pliku
void readData(std::string filename, upcxx::global_ptr<double> dataVector)
{
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
                        upcxx::rput(atof(token.c_str()), dataVector + count*7).wait();
                        break;
                    case 2:
                        upcxx::rput(atof(token.c_str()), dataVector + count*7 + 1).wait();
                        break;
                    case 3:
                        upcxx::rput(atof(token.c_str()), dataVector + count*7 + 2).wait();
                        break;
                    case 4:
                        upcxx::rput(atof(token.c_str()), dataVector + count*7 + 3).wait();
                        break;
                    case 5:
                        upcxx::rput(atof(token.c_str()), dataVector + count*7 + 4).wait();
                        break;
                    case 6:
                        upcxx::rput(atof(token.c_str()), dataVector + count*7 + 5).wait();
                        break;
                    }
                    line.erase(0, pos + delimiter.length());
                    datapos++;
                }
                upcxx::rput(atof(line.c_str()), dataVector + count*7 + 6).wait();
            }
            ++count;
        }
        datafile.close();
    }
}

//Inicjalizacja wektorow lokalnych danych na podstawie dataVector
void initializeVectors(int totalObjectCount, upcxx::global_ptr<double> dataVector, double* xPositionVector, double* yPositionVector, double* zPositionVector, double* xVelocityVector, double* yVelocityVector, double* zVelocityVector, double* massVector)
{
    for(int i = 0; i < totalObjectCount; i++)
    {
        xPositionVector[i] = upcxx::rget(dataVector + i*PROPERTIES).wait();
        yPositionVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 1).wait();
        zPositionVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 2).wait();
        xVelocityVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 3).wait();
        yVelocityVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 4).wait();
        zVelocityVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 5).wait();
        massVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 6).wait();
    }
}

//Aktualizacja dataVector na podstawie wektorow lokalnych
void updateDataVector(int i, upcxx::global_ptr<double> dataVector, double* xPositionVector, double* yPositionVector, double* zPositionVector)
{
    upcxx::rput(xPositionVector[i], dataVector + i*PROPERTIES).wait();
    upcxx::rput(yPositionVector[i], dataVector + i*PROPERTIES + 1).wait();
    upcxx::rput(zPositionVector[i], dataVector + i*PROPERTIES + 2).wait();
}

//Aktualizacja wektorow lokalnych z dataVector
void updatePositionVectors(int ownObjectStart, int ownObjectEnd, int totalObjectCount, upcxx::global_ptr<double> dataVector, double* xPositionVector, double* yPositionVector, double* zPositionVector)
{
    for(int i = 0; i < totalObjectCount; i++)
    {
        if(i < ownObjectStart || i > ownObjectEnd)
        {
            xPositionVector[i] = upcxx::rget(dataVector + i*PROPERTIES).wait();
            yPositionVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 1).wait();
            zPositionVector[i] = upcxx::rget(dataVector + i*PROPERTIES + 2).wait();
        }
    }
}

//Zapis danych z jednej iteracji do pliku wynikowego
void saveData(FILE* resultFile, double* xPositionVector, double* yPositionVector, double* zPositionVector, int totalObjectCount)
{
    for (int i = 0; i < totalObjectCount; i++)
    {
        fprintf(resultFile, "%d;%f;%f;%f\n", i, xPositionVector[i], yPositionVector[i], zPositionVector[i]);
    }
}