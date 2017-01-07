#include <QString>
#include <QtTest>

#include <RandLib.h>
#include "data.h"

class TestStableRand : public QObject
{
    Q_OBJECT

public:
    TestStableRand();

private Q_SLOTS:
    void testPdf1();
    void testPdf2();
    void testCdf();
};

TestStableRand::TestStableRand()
{
}

/**
 * @brief TestStableRand::testPdf
 * This test checks if implemented pdf of stable distribution
 * deviates from table values in major points less than on 0.1%.
 * Otherwise, the erreur admissible is less than 0.0001
 * (due to accuracy of table values)
 */
void TestStableRand::testPdf1()
{
    /// Check for each exponent = 0.25, 0.5, ..., 1.75,
    /// each skewness = -1.0, -0.75, ..., 1.0,
    /// and x = 0.0, 0.1, ..., 5.0
    StableRand rv(1.0, -0.25, 1.0, 0.0);
    for (int k = 0; k != 7; ++k) {
        float alpha = (k + 1) * 0.25;
        for (int i = 0; i != 9; ++i) {
            float beta = Pdf_data[k][i];
            rv.SetParameters(alpha, beta);
            for (int j = 0; j != 51; ++j) {
                float x = Pdf_data[k][j * 10 + 9];
                float x1_unrounded = rv.f(x);
                float x2 = Pdf_data[k][(j + 1) * 10 + i];
                float x1 = std::round(x1_unrounded * 1e4) / 1e4;
                bool isOk = false;
                if (x2 != 0.0) {
                    double fx1 = std::fabs(x1), fx2 = std::fabs(x2), fx1_unrounded = std::fabs(x1_unrounded);
                    double error1 = (std::fabs(x1 - x2) / std::max(fx1, fx2));
                    double error2 = (std::fabs(x1_unrounded - x2) / std::max(fx1_unrounded, fx2));
                    double error = std::min(error1, error2);
                    isOk = error < 0.001;
                }
                else {
                    isOk = (x1 == 0.0);
                }
                if (!isOk)
                    isOk = std::fabs(x1_unrounded - x2) < 1e-4;
                QVERIFY2(isOk, qPrintable("x = " + QString::number(x) + ", f(x) = " + QString::number(x1_unrounded) + " but should be: " + QString::number(x2) +
                                          " (alpha = " + QString::number(alpha) + " and beta = " + QString::number(beta) + ")"));
            }
        }
    }

}

/**
 * @brief TestStableRand::testPdf2
 * This test checks if implemented pdf of symmetric stable distribution
 * deviates from table values in major points less than on 0.1%.
 * Otherwise, the erreur admissible is less than 0.0001
 * (due to accuracy of table values)
 */
void TestStableRand::testPdf2()
{
    /// Check for each exponent = 1.0, 1.1, ..., 2.0,
    /// skewness is 0
    /// and x = 0.0, 0.1, ..., 1.0, 1.05, ..., 3.0
    StableRand rv(2.0, 0.0, 1.0, 0.0);
    for (int i = 0; i != 11; ++i) {
        float alpha = Pdf_data2[i];
        double alphaAdj = std::pow(alpha, 1.0 / alpha);
        rv.SetParameters(alpha, 0.0);
        for (int j = 0; j != 51; ++j) {
            float x = Pdf_data2[j * 12 + 11];
            float xAdj = alphaAdj * x;
            float y = rv.f(xAdj);
            float x1_unrounded = alphaAdj * y;
            float x2 = Pdf_data2[(j + 1) * 12 + i];
            float x1 = std::round(x1_unrounded * 1e4) / 1e4;
            bool isOk = false;
            if (x2 != 0.0) {
                double fx1 = std::fabs(x1), fx2 = std::fabs(x2), fx1_unrounded = std::fabs(x1_unrounded);
                double error1 = (std::fabs(x1 - x2) / std::max(fx1, fx2));
                double error2 = (std::fabs(x1_unrounded - x2) / std::max(fx1_unrounded, fx2));
                double error = std::min(error1, error2);
                isOk = error < 0.001;
            }
            else {
                isOk = (x1 == 0.0);
            }
            if (!isOk)
                isOk = std::fabs(x1_unrounded - x2) < 1e-4;
            QVERIFY2(isOk, qPrintable("x = " + QString::number(x) + ", f(x) = " + QString::number(y) + " but should be: " + QString::number(x2 / alphaAdj)
                                      + " and hf(hx) = " + QString::number(x1_unrounded) + " but should be: " + QString::number(x2)
                                      + " (alpha = " + QString::number(alpha) + " and beta = " + QString::number(0) + ")"));
        }
    }
}

void TestStableRand::testCdf()
{
    QVERIFY2(true, "f");
}

QTEST_APPLESS_MAIN(TestStableRand)

#include "tst_stable.moc"
