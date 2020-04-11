#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "distributionFunctions.h"

using namespace std;

int main() {
    ofstream file;

    //Time
    time_t now = time(0);
    char* dt = ctime(&now);
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);

    file.open ("../ANOVA-Results.txt");
    file << "===================================================================================" << endl;
    file << "  ANOVA i KONTRASTI                                       "<<dt;
    file << "===================================================================================" << endl;
    cout << "===================================================================================" << endl;
    cout << "  ANOVA i KONTRASTI                                       "<<dt;
    cout << "===================================================================================" << endl;
    int numberOfAlternatives = 0;
    int numberOfMeasurements = 0;
    double globalMean=0;
    double df1,df2;
    double significanceLevel = 0;


    //At least one alternative with two measurements
    while (numberOfAlternatives <= 1 || numberOfMeasurements <= 1) {
        cout << "Unesite broj alternativa: ";
        cin >> numberOfAlternatives;
        cout << "Unesite broj mjerenja za alternative: ";
        cin >> numberOfMeasurements;
        if(numberOfAlternatives <= 1 || numberOfMeasurements <= 1)
            cout<<"GREŠKA: Morate imati bar dvije alternative i bar dva mjerenja. Unesite ponovo."<<endl;
    }
    cout << "Unesite nivo povjerenja: ";
    cin >> significanceLevel;
    cout<<endl;
    //Create a mateix for storing data
    double measurementMatrix[numberOfAlternatives][numberOfMeasurements];
    df1 = numberOfAlternatives-1;
    df2 = numberOfAlternatives*(numberOfMeasurements-1);
    //Fill the matrix with data
    for(short i=0;i<numberOfAlternatives;i++)
    {
        cout<<"\rAlternativa "<<i+1<<" :"<<endl;
        cout<<"\r=============================="<<endl;
        for(short j=0;j<numberOfMeasurements;j++)
        {
            cout << "\r" << j + 1 << ". mjerenje: ";
            cin >> measurementMatrix[i][j];
        }
    }
    cout<<endl;
    cout<<"Generisanje rezultata može potrajati i do 30 sekundi... "<<endl;
    //Calculate mean of every alternative
    double columnMeans[numberOfAlternatives];
    for(short i=0;i<numberOfAlternatives;i++)
    {
        double mean=0,sum=0;
        for(short j=0;j<numberOfMeasurements;j++)
            sum+=measurementMatrix[i][j];
        mean = sum/(double)numberOfMeasurements;
        columnMeans[i]=mean;
    }
    //Calculate mean of all alternative means
    double sum=0;
    for(short i=0;i<numberOfAlternatives;i++)
        sum+=columnMeans[i];
    globalMean=sum/(double)numberOfAlternatives;

    //Calculate effects
    double effects[numberOfAlternatives];
    for(int i=0;i<numberOfAlternatives;i++)
        effects[i] = columnMeans[i] - globalMean;

    //Calculate SSA, SSE, SST
    double SSA=0,SSE=0, SST=0;
    //SSA
    sum=0;
    for(short i=0;i<numberOfAlternatives;i++)
        sum += pow(columnMeans[i]-globalMean,2);
    SSA = numberOfMeasurements * sum;
    //SSE
    for(short i=0;i<numberOfAlternatives;i++)
        for(short j=0;j<numberOfMeasurements;j++)
            SSE += pow(measurementMatrix[i][j]-columnMeans[i],2);
    //SST
    SST = SSA + SSE;

    //Calculate divisions
    double SSAdivSST = SSA/SST; //difference between systems
    double SSEdivSST = SSE/SST; //measurement errors

    //Variances
    double aVariance = SSA / (numberOfAlternatives-1);
    double eVariance = SSE / (numberOfAlternatives*(numberOfMeasurements-1));

    //Calculate F-value
    double F_value = aVariance/eVariance;
    double F_valueTable=0;


    file<<endl;
    file<<"==================================================================================="<<endl;
    file<<"|                                     ANOVA                                       |" << endl;
    file<<"==================================================================================="<<endl<<endl;

    //Print calculated table
    short cellWidth = (numberOfAlternatives+2)*10;

    for(short i=0;i<cellWidth;i++) file<<"-"; file<<endl;
    file<<"          ALTERNATIVE"<<endl;
    //for(short i=0;i<cellWidth;i++) cout<<"-"; cout<<endl;
    file<<"RB.MJER.  ";
    for(short i=0;i<numberOfAlternatives;i++) file<<i+1<<".        ";
    file<<"UKUP.SR.V.";
    file<<endl;
    for(short i=0;i<cellWidth;i++) file<<"-"; file<<endl;
    for(int i=0;i<numberOfMeasurements;i++)
    {
        file<<left<<setw(10)<<i+1;
        for(int j=0;j<numberOfAlternatives;j++)
            file<<left<<setw(10)<<fixed<<setprecision(5)<<measurementMatrix[j][i];
        file<<endl;
    }

    file<<"S.V.KOL.  ";
    for(short i=0;i<numberOfAlternatives;i++)
        file<<left<<setw(10)<<fixed<<setprecision(5)<<columnMeans[i];
    file<<left<<setw(10)<<fixed<<setprecision(5)<<globalMean;
    file<<endl;

    file<<"EFEKTI    ";
    for(short i=0;i<numberOfAlternatives;i++)
        file<<left<<setw(10)<<fixed<<setprecision(5)<<effects[i];
    file<<endl;
    for(short i=0;i<cellWidth;i++) file<<"-"; file<<endl;

    file<<"SSA       "<<left<<setw(10)<<fixed<<setprecision(5)<<SSA<<endl;
    file<<"SSE       "<<left<<setw(10)<<fixed<<setprecision(5)<<SSE<<endl;
    file<<"SST       "<<left<<setw(10)<<fixed<<setprecision(5)<<SST<<endl;

    file<<"SSA/SST   "<<left<<setw(10)<<fixed<<setprecision(5)<<SSAdivSST<<endl;
    file<<"SSE/SST   "<<left<<setw(10)<<fixed<<setprecision(5)<<SSEdivSST<<endl;

    file<<"SSE var   "<<left<<setw(10)<<fixed<<setprecision(5)<<eVariance<<endl;
    file<<"SSA var   "<<left<<setw(10)<<fixed<<setprecision(5)<<aVariance<<endl;

    file<<"F[izrač.] "<<left<<setw(10)<<fixed<<setprecision(5)<<F_value<<endl;
    F_valueTable = fDistributionTableValue(df1,df2,significanceLevel);
    file<<"F[tabel.] "<<left<<setw(10)<<fixed<<setprecision(3)<<F_valueTable<<endl;
    for(short i=0;i<cellWidth;i++) file<<"-"; file<<endl;

    file<<"Na osnovu dobijenih vrijednosti može se reći da je "<<SSAdivSST*100<<"% ukupne varijacije zbog"<<endl;
    file<<"razlika između sistema, dok je "<<SSEdivSST*100<<"% ukupne varijacije zbog grešaka u mjerenju."<<endl;
    if(F_value < F_valueTable)
    {
        file<<"Pošto je izračunara F-vrijednost manja od F-vrijednosti dobijene iz tabele"<<endl;
        file<<"može se zaključiti da sa nivoom povjerenja od "<<setprecision(1)<<(1-significanceLevel)*100<<"% ne postoji"<<endl;
        file<<"statistički značajna razlika između sistema."<<endl;
    }
    else
    {
        file<<"Pošto je izračunara F-vrijednost veća od F-vrijednosti dobijene iz tabele"<<endl;
        file<<"može se zaključiti da sa nivoom povjerenja od "<<setprecision(1)<<(1-significanceLevel)*100<<"% postoji"<<endl;
        file<<"statistički značajna razlika između sistema."<<endl;
    }
    for(short i=0;i<cellWidth;i++) file<<"-"; file<<endl;

    file<<endl;
    file<<"==================================================================================="<<endl;
    file<<"|                                     KONTRASTI                                   |" << endl;
    file<<"==================================================================================="<<endl<<endl;


    for(int i=0;i<numberOfAlternatives;i++)
    {
        for (int j=1+i;j<numberOfAlternatives;j++)
        {
            double contrast = effects[i]-effects[j];
            double var = (sqrt(eVariance)) * (sqrt(2/(double)(numberOfAlternatives*numberOfMeasurements)));
            double tValue = tDistributionTableValue(df2, significanceLevel);
            //Calculate intervals
            double c1 = contrast + tValue * var;
            double c2 = contrast - tValue * var;


            file<<"----------------------------------------------------------------------------------"<<endl;
            file<< "POREĐENJE ALTERNATIVA "<<i+1<<"-"<<j+1<<endl;
            file<<"----------------------------------------------------------------------------------"<<endl;
            file<<"Kontrast            "<<left<<setw(20)<<fixed<<setprecision(10)<<contrast<<endl;
            file<<"Varijansa           "<<left<<setw(20)<<fixed<<setprecision(10)<<var<<endl;
            file<<"t-vrijednost        "<<left<<setw(20)<<fixed<<setprecision(10)<<tValue<<endl;
            file<<"Interval povjerenja "<<left<<fixed<<setprecision(7)<<"["<<c2<<","<<c1<<"]"<<endl;
            if(c1>0 && c2<0)
            {
                file<<"Pošto ovaj interval povjerenja uključuje nulu, zaključuje se da ne postoji"<<endl;
                file<<"statistički značajna razlika između alternativa."<<endl;
            }
            else
            {
                file<<"Pošto ovaj interval povjerenja ne uključuje nulu, može se zaključiti da se"<<endl;
                file<<"ne može reći da ne postoji značajna razlika između ove dvije alternative."<<endl;
            }
        }
    }

    file<<"==================================================================================="<<endl;
    file<<"|                              ©Želimir Maletić, 2020                             |" << endl;
    file<<"|                              mzeljo417@outlook.com                              |" << endl;
    file<<"==================================================================================="<<endl;

    cout<<"Analiza završena!"<<endl;
    cout<<"Rezultati su generisani i zapisani u fajl ANOVA-Results.txt !"<<endl;
    file.close();
    return 0;
}
