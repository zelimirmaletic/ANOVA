#ifndef ANOVA_DISTRIBUTIONFUNCTIONS_H
#define ANOVA_DISTRIBUTIONFUNCTIONS_H
double fDistributionPDF(double x, double df1, double df2);
double tDistributionPDF(double x, double df);
double beta(double x, double y);
double  simpsonsFormula(double a, double b, int m);
double fDistributionTableValue(double degreesOfFreedom1, double degreesOfFreedom2, double significanceLevel);
double tDistributionTableValue(double degreesOfFreedom, double significanceLevel);
#endif //ANOVA_DISTRIBUTIONFUNCTIONS_H
