#include "RandMath.h"

namespace RandMath
{

bool areClose(double a, double b, double eps)
{
    if (a == b)
        return true;
    double fa = std::fabs(a);
    double fb = std::fabs(b);
    return (std::fabs(b - a) < eps * std::max(fa, fb));
}

int sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double atan(double x)
{
    /// For small absolute values we use standard technique
    /// Otherwise we use relation
    /// atan(x) = +/-π/2 - atan(1/x)
    /// to avoid numeric problems
    if (x == 0.0)
        return 0.0;
    if (x > 1.0)
        return M_PI_2 - std::atan(1.0 / x);
    return (x < -1.0) ? -M_PI_2 - std::atan(1.0 / x) : std::atan(x);
}

double log1pexp(double x)
{
    if (x < 20.0)
        return std::log1p(std::exp(x));
    return (x < 35.0) ? x + std::exp(-x) : x;
}

double log1mexp(double x)
{
    return (x < -M_LN2) ? std::log1p(-std::exp(x)) : std::log(-std::expm1(x));
}

double logexpm1(double x)
{
    if (x < 20.0)
        return std::log(std::expm1(x));
    return (x < 35.0) ? x - std::exp(-x) : x;
}

double log2mexp(double x)
{
    return std::log1p(-std::expm1(x));
}

/**
 * @brief FACTORIAL_TABLESIZE maximum value for input parameter to use table method
 */
constexpr int FACTORIAL_TABLESIZE = 255;

/**
 * @brief FACTORIAL_TABLE (n * 10)! for n from 0 to 25
 */
constexpr long double FACTORIAL_TABLE[] =
{
    1.l,
    3628800.l,
    2432902008176640000.l,
    265252859812191058636308480000000.l,
    815915283247897734345611269596115894272000000000.l,
    30414093201713378043612608166064768844377641568960512000000000000.l,
    8320987112741390144276341183223364380754172606361245952449277696409600000000000000.l,
    11978571669969891796072783721689098736458938142546425857555362864628009582789845319680000000000000000.l,
    71569457046263802294811533723186532165584657342365752577109445058227039255480148842668944867280814080000000000000000000.l,
    1485715964481761497309522733620825737885569961284688766942216863704985393094065876545992131370884059645617234469978112000000000000000000000.l,
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000.l,
    15882455415227429404253703127090772871724410234473563207581748318444567162948183030959960131517678520479243672638179990208521148623422266876757623911219200000000000000000000000000.l,
    6689502913449127057588118054090372586752746333138029810295671352301633557244962989366874165271984981308157637893214090552534408589408121859898481114389650005964960521256960000000000000000000000000000.l,
    6466855489220473672507304395536485253155359447828049608975952322944781961185526165512707047229268452925683969240398027149120740074042105844737747799459310029635780991774612983803150965145600000000000000000000000000000000.l,
    13462012475717524605876073858941615558355851148193967190051391468057460367090535696797920946629681836680869097041958983702264048370902871114013579941370766400374327741701139895604871545254810788060989321379840000000000000000000000000000000000.l,
    57133839564458545904789328652610540031895535786011264182548375833179829124845398393126574488675311145377107878746854204162666250198684504466355949195922066574942592095735778929325357290444962472405416790722118445437122269675520000000000000000000000000000000000000.l,
    471472363599206132240694321176194377951192623045460204976904578317542573467421580346978030238114995699562728104819596262106947389303901748942909887857509625114880781313585012959529941660203611234871833992565791817698209861793313332044734813700096000000000000000000000000000000000000000.l,
    7257415615307998967396728211129263114716991681296451376543577798900561843401706157852350749242617459511490991237838520776666022565442753025328900773207510902400430280058295603966612599658257104398558294257568966313439612262571094946806711205568880457193340212661452800000000000000000000000000000000000000000.l,
    200896062499134299656951336898466838917540340798867777940435335160044860953395980941180138112097309735631594101037399609671032132186331495273609598531966730972945653558819806475064353856858157445040809209560358463319644664891114256430017824141796753818192338642302693327818731986039603200000000000000000000000000000000000000000000.l,
    9680322675255249156123346514615331205418161260462873360750859919944104623425228207640470674933540169424682360525991982916161596983449594045525553704253602287443197783274656957056546338783001340434094795097553229620273057440272298773179365935914105128629426348958748638226084106818484328004851174161755668480000000000000000000000000000000000000000000000.l,
    788657867364790503552363213932185062295135977687173263294742533244359449963403342920304284011984623904177212138919638830257642790242637105061926624952829931113462857270763317237396988943922445621451664240254033291864131227428294853277524242407573903240321257405579568660226031904170324062351700858796178922222789623703897374720000000000000000000000000000000000000000000000000.l,
    105823620292236563784274284243348353057589905787169019562352737522144487532400210147849369011714673954768265316577892528273760626189481169051055226066650741189573897273684791411180134039439160066561895838501000817711682625725670477616267598661259194975646029749546282594356217374097544153589482020891750774735012558313460846824864172030239122128896000000000000000000000000000000000000000000000000000.l,
    22838603359146414573972658651153337270429730715462287017736347161260276926030248458777765497919211029457065581960747795750095505232241970499561769723020565876672261660609763234049775547325430135571331468257475537994508495233770658945310210552725163342784668756149049213658078338458534285571551800849578848226429898670032945513859929938621783523490272646966918544936140800000000000000000000000000000000000000000000000000000.l,
    7758587304686725201813174298892781442413952130995533365303964524344944412641389739603152000644515957408814002319492032321234250506968028455594445689972313374305301019340949789291189972149450405025159624155827152329676580440959428615802893638146558163235483142136540783687811997927615346859658417205832954125915861983307177232587595821512723429698627780530255874167602077755356592824804966400000000000000000000000000000000000000000000000000000000.l,
    4067885363647058120493575921486885310172051259182827146069755969081486918925585104009100729728348522923820890245870098659147156051905732563147381599098459244752463027688115705371704628286326621238456543307267608612545168337779669138759451760395968217423617954330737034164596496963986817722252221059768080852489940995605579171999666916004042965293896799800598079985264195119506681577622056215044851618236292136960000000000000000000000000000000000000000000000000000000000.l,
    3232856260909107732320814552024368470994843717673780666747942427112823747555111209488817915371028199450928507353189432926730931712808990822791030279071281921676527240189264733218041186261006832925365133678939089569935713530175040513178760077247933065402339006164825552248819436572586057399222641254832982204849137721776650641276858807153128978777672951913990844377478702589172973255150283241787320658188482062478582659808848825548800000000000000000000000000000000000000000000000000000000000000.l,
};

/**
 * @brief factorialForSmallValue
 * Get n! using product method with table
 * @param n non-negative integer number
 * @return n!
 */
long double factorialForSmallValue(int n)
{
    int residue = n % 10;
    if (residue <= 5)
    {
        /// go up
        int nPrev = n - residue;
        long double fact = FACTORIAL_TABLE[nPrev / 10];
        for (int i = 1; i <= residue; ++i)
            fact *= nPrev + i;
        return fact;
    }

    /// go  down
    int nNext = n - residue + 10;
    double denominator = 1;
    for (int i = 0; i < 10 - residue; ++i)
        denominator *= nNext - i;
    return FACTORIAL_TABLE[nNext / 10] / denominator;
}

long double factorial(double n)
{
    if (n < 0)
        return 0.0;
    return (n > FACTORIAL_TABLESIZE) ? std::tgamma(n + 1) : factorialForSmallValue(n);
}

long double doubleFactorial(int n)
{
    long double n_fact = factorial(n);
    if (n & 1) {
        n <<= 1;
        return factorial(n + 1) / (n * n_fact);
    }
    return (1 << n) * n_fact;
}

long double binomialCoef(int n, int k)
{
    if (k > n)
        return 0;
    if (n > 20) /// to avoid overflow
        return 1.0 / ((n + 1) * beta(n - k + 1, k + 1));
    long double n_fact = factorial(n);
    long double k_fact = factorial(k);
    long double n_k_fact;
    if (k == n - k)
        n_k_fact = k_fact;
    else
        n_k_fact = factorial(n - k);
    return n_fact / (k_fact * n_k_fact);
}

double erfinvChebyshevSeries(double x, long double t, const long double *array, int size)
{
    /// We approximate inverse erf via Chebyshev polynomials
    long double Tn = t, Tnm1 = 1, sum = 0.0;
    for (int i = 1; i != size; ++i) {
        sum += array[i] * Tn;
        /// calculate next Chebyshev polynomial
        double temp = Tn;
        Tn *= 2 * t;
        Tn -= Tnm1;
        Tnm1 = temp;
    }
    return x * (array[0] + sum);
}

double erfinvAux1(double beta)
{
    /// |1 - p| < 5e-16
    static constexpr long double D5 = -9.199992358830151031278420l, D6 = 2.794990820124599493768426l;
    static constexpr int muSize = 26;
    static constexpr long double mu[muSize] = {.9885750640661893136460358l, .0108577051845994776160281l, -.0017511651027627952594825l, .0000211969932065633437984l,
                                               .0000156648714042435087911l, -.05190416869103124261e-5l, -.00371357897426717780e-5l, .00012174308662357429e-5l,
                                              -.00001768115526613442e-5l, -.119372182556161e-10l, .003802505358299e-10l, -.000660188322362e-10l, -.000087917055170e-10l,
                                              -.3506869329e-15l, -.0697221497e-15l, -.0109567941e-15l, -.0011536390e-15l, -.0000263938e-15l, .05341e-20l, -.22610e-20l,
                                               .09552e-20l, -.05250e-20l, .02487e-20l, -.01134e-20l, .00420e-20l};
    return erfinvChebyshevSeries(beta, D5 / std::sqrt(beta) + D6, mu, muSize);
}

double erfinvAux2(double beta)
{
    /// 5e-16 < |1 - p| < 0.025
    static constexpr long double D3 = -0.5594576313298323225436913l, D4 = 2.287915716263357638965891l;
    static constexpr int deltaSize = 38;
    static constexpr long double delta[deltaSize] = {.9566797090204925274526373l, -.0231070043090649036999908l, -.0043742360975084077333218l, -.0005765034226511854809364l,
                                                    -.0000109610223070923931242l, .0000251085470246442787982l, .0000105623360679477511955l, .27544123300306391503e-5l,
                                                     .04324844983283380689e-5l, -.00205303366552086916e-5l, -.00438915366654316784e-5l, -.00176840095080881795e-5l,
                                                    -.00039912890280463420e-5l, -.00001869324124559212e-5l, .00002729227396746077e-5l, .00001328172131565497e-5l,
                                                     .318342484482286e-10l, .016700607751926e-10l, -.020364649611537e-10l, -.009648468127965e-10l, -.002195672778128e-10l,
                                                    -.000095689813014e-10l, .000137032572230e-10l, .000062538505417e-10l, .000014584615266e-10l, .1078123993e-15l,
                                                    -.0709229988e-15l, -.0391411775e-15l, -.0111659209e-15l, -.0015770366e-15l, .0002853149e-15l, .0002716662e-15l,
                                                     .0000176835e-15l, .09828e-20l, .20464e-20l, .08020e-20l, .01650e-20l};
    return erfinvChebyshevSeries(beta, D3 * beta + D4, delta, deltaSize);
}

double erfinvAux3(double beta)
{
    /// 0.8 < p < 0.975
    static constexpr long double D1 = -1.548813042373261659512742l, D2 = 2.565490123147816151928163l;
    static constexpr int lambdaSize = 27;
    static constexpr long double lambda[lambdaSize] = {.9121588034175537733059200l, -.0162662818676636958546661l, .0004335564729494453650589l, .0002144385700744592065205l,
                                                       .26257510757648130176e-5l, -.30210910501037969912e-5l, -.00124060618367572157e-5l, .00624066092999917380e-5l,
                                                      -.00005401247900957858e-5l, -.00014232078975315910e-5l, .343840281955305e-10l, .335848703900138e-10l, -.014584288516512e-10l,
                                                      -.008102174258833e-10l, .000525324085874e-10l, .000197115408612e-10l, -.000017494333828e-10l, -.4800596619e-15l, .0557302987e-15l,
                                                       .0116326054e-15l, -.0017262489e-15l, -.0002784973e-15l, .0000524481e-15l, .65270e-20l, -.15707e-20l, -.01475e-20l, .00450e-20l};
    return erfinvChebyshevSeries(beta, D1 * beta + D2, lambda, lambdaSize);
}

double erfinvAux4(double p)
{
    /// 0 < p < 0.8
    static constexpr int xiSize = 39;
    static constexpr long double xi[xiSize] = {.9928853766189408231495800l, .1204675161431044864647846l, .0160781993420999447267039l, .0026867044371623158279591l,
                                               .0004996347302357262947170l, .0000988982185991204409911l, .0000203918127639944337340l, .43272716177354218758e-5l,
                                               .09380814128593406758e-5l, .02067347208683427411e-5l, .00461596991054300078e-5l, .00104166797027146217e-5l,
                                               .00023715009995921222e-5l, .00005439284068471390e-5l, .00001255489864097987e-5l, .291381803663201e-10l,
                                               .067949421808797e-10l, .015912343331569e-10l, .003740250585245e-10l, .000882087762421e-10l, .000208650897725e-10l,
                                               .000049488041039e-10l, .000011766394740e-10l, .2803855725e-15l, .0669506638e-15l, .0160165495e-15l, .0038382583e-15l,
                                               .0009212851e-15l, .0002214615e-15l, .0000533091e-15l, .0000128488e-15l, .31006e-20l, .07491e-20l, .01812e-20l,
                                               .00439e-20l, .00106e-20l, .00026e-20l, .00006e-20l, .00002e-20l};
    return erfinvChebyshevSeries(p, p * p / 0.32 - 1.0, xi, xiSize);
}

double erfinv(double p)
{
    /// Consider special cases
    if (p < 0.0)
        return -erfinv(-p);
    if (p > 1.0)
        return NAN;
    if (p == 1.0)
        return INFINITY;
    if (p == 0.0)
        return 0.0;
    if (p < 0.8)
        return erfinvAux4(p);
    /// Handle tails
    double beta = std::sqrt(-std::log1p(-p * p));
    if (p < 0.9975)
        return erfinvAux3(beta);
    return (1.0 - p < 5e-16) ? erfinvAux1(beta) : erfinvAux2(beta);
}

double erfcinv(double p)
{
    /// Consider special cases
    if (p > 1.0)
        return -erfcinv(2.0 - p);
    if (p == 0.0)
        return INFINITY;
    if (p == 1.0)
        return 0.0;
    if (p > 0.2)
        return erfinvAux4(1.0 - p);
    double pSq = p * p, p2 = 2 * p;
    double beta = std::sqrt(-std::log(p2 - pSq));
    if (p > 0.0025)
        return erfinvAux3(beta);
    return (p > 5e-16) ? erfinvAux2(beta) : erfinvAux1(beta);
}

double linearInterpolation(double a, double b, double fa, double fb, double x)
{
    if (b == a)
        return fa;
    double fx = x - a;
    fx /= (b - a);
    fx *= (fb - fa);
    return fx + fa;
}

double harmonicNumber(double exponent, int number)
{
    if (number < 1)
        return 0;
    if (exponent == 1)
        return M_EULER + digamma(number + 1);
    if (exponent == 2)
        return M_PI_SQ / 6.0 - trigamma(number + 1);
    double res = 1.0;
    for (int i = 2; i <= number; ++i)
        res += std::pow(i, -exponent);
    return res;
}

double logModifiedBesselFirstKind(double x, double nu)
{
    // TODO: in c++17, use implemented series for too small (x ~< nu / 2) and too large x > ?
    // otherwise, it is better to use log(besseli(x, n))
    if (x < 0) {
        double roundNu = std::round(nu);
        bool nuIsInt = areClose(nu, roundNu);
        if (nuIsInt) {
            int nuInt = roundNu;
            return (nuInt % 2) ? NAN : logModifiedBesselFirstKind(-x, nu);
        }
        return -INFINITY;
    }

    if (x == 0) {
        if (nu == 0)
            return 0.0;
        double roundNu = std::round(nu);
        bool nuIsInt = areClose(nu, roundNu);
        return (nu > 0 || nuIsInt) ? -INFINITY : INFINITY;
    }

    if (std::fabs(nu) == 0.5) {
        /// log(sinh(x)) or log(cosh(x))
        double y = 0.5 * std::log(M_2_PI / x);
        y -= M_LN2;
        y += x;
        y += (nu > 0) ? RandMath::log1pexp(-2 * x) : RandMath::log1mexp(-2 * x);
        return y;
    }

    /// small x
    if (x < 2 * nu)
    {
        double addon = 1.0, sum = 0.0;
        int i = 1;
        double logHalfx = std::log(0.5 * x);
        do {
            addon = 2 * i * logHalfx;
            addon -= std::lgamma(nu + i + 1);
            addon = std::exp(addon - lfact(i));
            sum += addon;
            ++i;
        } while (std::fabs(addon) > MIN_POSITIVE * std::fabs(sum));
        double y = std::log(sum + 1.0 / std::tgamma(nu + 1));
        y += nu * logHalfx;
        return y;
    }

    double addon = 1.0;
    double sum = 1.0;
    double denominator = 0.125 / x;
    double n2Sq = 4.0 * nu * nu;
    double i = 1.0, j = 1.0;
    do {
        double numerator = j * j - n2Sq;
        double frac = numerator * denominator / i;
        if (frac > 1)
            break;
        addon *= frac;
        sum += addon;
        ++i;
        j += 2;
    } while (std::fabs(addon) > MIN_POSITIVE * std::fabs(sum));
    double y = x;
    y -= 0.5 * std::log(2 * M_PI * x);
    y += std::log(sum);
    return y;
}

double logModifiedBesselSecondKind(double x, double nu)
{
    if (nu == 0.5)
        return 0.5 * std::log(M_PI_2 / x) - x;

    double y = std::exp(logModifiedBesselFirstKind(x, -nu)); // besseli
    y -= std::exp(logModifiedBesselFirstKind(x, nu)); // besseli
    y /= std::sin(nu * M_PI);
    y *= M_PI_2;
    return std::log(y);
}

double zetaRiemann(double s)
{
    if (s == 1)
        return INFINITY;
    if (s < 1)
        return NAN;
    int N = 100;
    double y = harmonicNumber(s, N);
    double NS = std::pow(N, -s);
    return y + N * NS / (s - 1) + 0.5 * NS;
}

/**
 * @brief WLambert
 * @param x
 * @param w0
 * @param epsilon
 * @return
 */
double WLambert(double x, double w0, double epsilon)
{
    double w = w0;
    double step = 0;
    do {
        double ew = std::exp(w);
        double wew = w * ew;
        double numerator1 = wew - x;
        double wp1 = w + 1;
        double denominator1 = ew * wp1;
        double numerator2 = (w + 2) * numerator1;
        double denominator2 = 2 * wp1;
        step = numerator2 / denominator2;
        step = numerator1 / (denominator1 - step);
        w -= step;
    } while (std::fabs(step) > epsilon);
    return w;
}

double W0Lambert(double x, double epsilon)
{
    double w = 0;
    if (x < -M_1_E)
        return NAN;
    if (x > 10) {
        double logX = std::log(x);
        double loglogX = std::log(logX);
        w = logX - loglogX;
    }
    return WLambert(x, w, epsilon);
}

double Wm1Lambert(double x, double epsilon)
{
    double w = -2;
    if (x < -M_1_E || x > 0)
        return NAN;
    if (x > -0.1) {
        double logmX = std::log(-x);
        double logmlogmX = std::log(-logmX);
        w = logmX - logmlogmX;
    }
    return WLambert(x, w, epsilon);
}

/**
 * @brief MarcumPSeries
 * @param mu
 * @param x
 * @param y
 * @return series expansion for Marcum-P function
 */
double MarcumPSeries(double mu, double x, double y)
{
    static constexpr double ln2piEps = -35.0; /// ~log(2πε) for ε = 1e-16
    double lgammamu = std::lgamma(mu);
    double C = lgammamu - ln2piEps + mu;
    double logX = std::log(x), logY = std::log(y);

    /// solving equation f(n) = 0
    /// to specify first negleted term
    double root = std::max(0.5 * (mu * mu + 4 * x * y - mu), 1.0);
    double logXY = logX + logY;
    if (!RandMath::findRoot([C, mu, logXY] (double n)
    {
        double npmu = n + mu;
        double logn = std::log(n), lognpmu = std::log(npmu);
        double first = logn - 2 - logXY;
        first *= n;
        first += npmu * lognpmu;
        first -= C;
        double second = logn + lognpmu - logXY;
        double third = 2.0 / n;
        return DoubleTriplet(first, second, third);
    }, root))
        return NAN; /// unexpected return

    /// series expansion
    double sum = 0.0;
    int n0 = std::max(std::ceil(root), 5.0); /// sanity check
    double P = pgamma(mu + n0, y, logY);
    for (int n = n0; n > 0; --n) {
        double term = n * logX - x;
        term = std::exp(term - lfact(n)) * P;
        sum += term;
        // check if diffP can be calculated via recursion
        double diffP = (mu + n - 1) * logY - y - std::lgamma(mu + n);
        P += std::exp(diffP);
    }
    sum += std::exp(-x) * P;
    return sum;
}

/**
 * @brief MarcumPAsymptoticForLargeXY
 * @param mu
 * @param x
 * @param y
 * @param sqrtX
 * @param sqrtY
 * @return asymptotic expansion for Marcum-P function for large x*y
 */
double MarcumPAsymptoticForLargeXY(double mu, double x, double y, double sqrtX, double sqrtY)
{
    double xi = 2 * sqrtX * sqrtY;
    double sigma = (y + x) / xi - 1;
    double sigmaXi = sigma * xi;
    double rho = sqrtY / sqrtX;
    double aux = std::erfc(sqrtX - sqrtY);
    double Phi = (sigma == 0) ? 0.0 : sign(x - y) * std::sqrt(M_PI / sigma) * aux;
    double Psi0 = aux / std::sqrt(rho);
    double logXi = M_LN2 + 0.5 * std::log(x * y);
    int n0 = std::max(std::ceil(sigmaXi), 7.0); /// sanity check
    double sum = 0.0;
    double A1 = 1.0, A2 = 1.0;
    for (int n = 1; n <= n0; ++n) {
        /// change φ
        double nmHalf = n - 0.5;
        Phi *= -sigma;
        double temp = -sigmaXi - nmHalf * logXi;
        Phi += std::exp(temp);
        Phi /= nmHalf;
        /// calculate A(μ-1) and A(μ)
        double coef1 = (nmHalf - mu) * (nmHalf + mu);
        double coef2 = (nmHalf - mu + 1) * (nmHalf + mu - 1);
        double denom = -2 * n;
        A1 *= coef1 / denom;
        A2 *= coef2 / denom;
        /// compute term ψ and add it to the sum
        double Psi = Phi * (A2 - A1 / rho);
        sum += (n & 1) ? Psi : -Psi;
    }
    sum /= M_SQRT2PI;
    sum += Psi0;
    sum *= 0.5 * std::pow(rho, mu);
    return sum;
}

/**
 * @brief MarcumPForMuLessThanOne
 * @param mu
 * @param x
 * @param y
 * @return
 */
double MarcumPForMuLessThanOne(double mu, double x, double y)
{
    // TODO: check Krishnamoorthy paper for alternative representation

    /// in this case we use numerical integration,
    /// however we have singularity point at 0,
    /// so we get rid of it by subtracting the function
    /// which has the same behaviour at this point
    double aux = x + mu * M_LN2 + std::lgamma(mu);
    double log2x = std::log(2 * x);
    double I = std::log(2 * y) * mu - aux;
    I = std::exp(I) / mu;

    double mum1 = mu - 1.0;
    I += RandMath::integral([x, log2x, mum1, aux] (double t)
    {
        if (t <= 0)
            return 0.0;
        /// Calculate log of leveling factor
        double log2T = std::log(2 * t);
        double exponent = mum1 * log2T;
        exponent -= aux;

        /// Calculate log(2f(t))
        double logBessel = RandMath::logModifiedBesselFirstKind(2 * std::sqrt(x * t), mum1);
        double z = mum1 * (log2T - log2x);
        double log2F = 0.5 * z - t - x + logBessel;

        /// Return difference f(t) - factor
        return std::exp(log2F) - 2 * std::exp(exponent);
    }, 0, y);
    return I;
}

double MarcumP(double mu, double x, double y)
{
    // TODO: hash values of sqrt(x), sqrt(y), log(x), log(y) if possible
    if (x < 0.0 || y <= 0.0)
        return 0.0;

    if (mu < 1.0)
        return MarcumPForMuLessThanOne(mu, x, y);

    if (x < 30)
        return MarcumPSeries(mu, x, y);

    double sqrtX = std::sqrt(x), sqrtY = std::sqrt(y);
    double xi = 2 * sqrtX * sqrtY;
    if (xi > 30 && mu * mu < 2 * xi)
        return MarcumPAsymptoticForLargeXY(mu, x, y, sqrtX, sqrtY);

    // TODO: implement the rest techniques

    double mum1 = mu - 1;
    double logX = std::log(x);
    return RandMath::integral([mum1, logX, x](double t){
        if (t < 0.0)
            return 0.0;
        if (t == 0.0)
            return (mum1 == 0) ? 0.5 * std::exp(-x) : 0.0;
        double y = RandMath::logModifiedBesselFirstKind(2 * std::sqrt(x * t), mum1);
        double z = 0.5 * mum1 * (std::log(t) - logX);
        double h = t + x;
        return std::exp(y + z - h);
    }, 0, y);
}

double MarcumQ(double mu, double x, double y)
{
    // TODO: implement and use, when mu + x > y
    return 1.0 - MarcumP(mu, x, y);
}

}

