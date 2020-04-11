#include "distributionFunctions.h"
#include <cmath>



double beta(double x, double y) {
    return ((tgamma(x)*tgamma(y))/tgamma(x+y));
}

double fDistributionPDF(double x, double df1, double df2)
{
    double a = (1/beta(df1/2,df2/2));
    double b = pow(df1/df2, df1/2);
    double c = pow(x, ((df1/2)-1));
    double d = pow(1+(df1/df2)*x, (-1*((df1+df2)/2)));
    return a*b*c*d;
}

double tDistributionPDF(double x, double df)
{
    double a = tgamma((df+1)/2);
    double b = sqrt(df*3.1415)*tgamma(df/2);
    double c = pow(1+((x*x)/df),(-1)*((df+1)/2));
    return (a/b)*c;
}

double  simpsonsFormulaF(double a, double b, int m, double df1, double df2)
{
        double h = (b-a)/(2*(double)m);
        double s1=0;
        double s2=0;
        double x=0;

        for(int k=1; k<m-1; k++)
        {
            x=a+h*(2*k);
            s1+= fDistributionPDF(x,df1,df2);
        }

        for(int k=1; k<m; k++)
        {
            x=a+h*(2*k-1);
            s2+= fDistributionPDF(x,df1,df2);
        }
        return (h*(fDistributionPDF(a,df1,df2)+fDistributionPDF(b,df1,df2)+2*s1+4*s2)/3);
}

double fDistributionTableValue(double degreesOfFreedom1, double degreesOfFreedom2, double significanceLevel)
{
    double tableValue = 1.6;
    double areaUnderCurve = 0.0;
    double precisionStep = 0.001;
    while(areaUnderCurve<1-significanceLevel)
    {
        areaUnderCurve = simpsonsFormulaF(0,tableValue,1000, degreesOfFreedom1, degreesOfFreedom2);
        tableValue+=precisionStep;
    }
    return tableValue;
}

double  simpsonsFormulaT(double a, double b, int m, double df)
{
    double h = (b-a)/(2*(double)m);
    double s1=0;
    double s2=0;
    double x=0;

    for(int k=1; k<m-1; k++)
    {
        x=a+h*(2*k);
        s1+= tDistributionPDF(x, df);
    }

    for(int k=1; k<m; k++)
    {
        x=a+h*(2*k-1);
        s2+= tDistributionPDF(x, df);
    }
    return (h*(tDistributionPDF(a, df)+tDistributionPDF(b, df)+2*s1+4*s2)/3);
}

double tDistributionTableValue(double degreesOfFreedom, double significanceLevel)
{
    double tableValue = 0.0;
    double areaUnderCurve = 0.0;
    double precisionStep = 0.001;
    double condition = 0.5-(significanceLevel);
    while(areaUnderCurve<condition)
    {
        areaUnderCurve = simpsonsFormulaT(0,tableValue,1000, degreesOfFreedom);
        tableValue+=precisionStep;
    }
    return tableValue;
}


