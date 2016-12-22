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

/**
 * @brief maxFactorialTableValue maximum value for input parameter to use table methods
 */
constexpr int maxFactorialTableValue = 255;

/**
 * @brief factorialTable (n * 10)! for n from 0 to 25
 */
constexpr long double factorialTable[] =
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
        long double fact = factorialTable[nPrev / 10];
        for (int i = 1; i <= residue; ++i)
            fact *= nPrev + i;
        return fact;
    }

    /// go  down
    int nNext = n - residue + 10;
    double denominator = 1;
    for (int i = 0; i < 10 - residue; ++i)
        denominator *= nNext - i;
    return factorialTable[nNext / 10] / denominator;
}

long double factorial(double n)
{
    if (n < 0)
        return 0.0;
    return (n > maxFactorialTableValue) ? std::tgamma(n + 1) : factorialForSmallValue(n);
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
    long double n_fact = factorial(n);
    long double k_fact = factorial(k);
    long double n_k_fact;
    if (k == n - k)
        n_k_fact = k_fact;
    else
        n_k_fact = factorial(n - k);
    return n_fact / (k_fact * n_k_fact);
}

double digamma(double x)
{
    /// Negative argument
    if (x < 0.0)
        return digamma(1.0 - x) + M_PI / std::tan(M_PI * (1.0 - x));

    /// Large argument
    if (x > 1000.0)
        return std::log(x) - 0.5 / x;

    /// Value at zero is undefined
    if (x == 0.0)
        return NAN; /// +/- INFINITY

    double y;
    if (x == 1.0)
        y = 0;
    else if (x == 0.5)
        y = 2 * M_LN2;
    else if (x == 0.25)
        y = -0.5 * M_PI - 3 * M_LN2;
    else {
        /// Use property: digamma(x) = digamma(x + 1) - 1 / x
        /// and integral representation
        y = integral([x] (double t)
        {
            if (t >= 1)
                return x;
            if (t <= 0)
                return 0.0;
            return (1.0 - std::pow(t, x)) / (1.0 - t);
        }, 0, 1) - 1.0 / x;
    }
    return y - M_EULER;
}

double trigamma(double x)
{
    /// Negative argument
    if (x < 0.0)
    {
        double z = M_PI / std::sin(M_PI * (1.0 - x));
        return z - digamma(1.0 - x);
    }

    /// Large argument
    if (x > 200.0)
        return (x + 0.5) / (x * x);

    /// Special values
    if (x == 0.0)
        return INFINITY;
    if (x == 0.25)
        return M_PI_SQ + 8 * M_CATALAN;
    if (x == 0.5)
        return 0.5 * M_PI_SQ;
    if (x == 1.0)
        return M_PI_SQ / 6.0;
    if (x == 1.5)
        return 0.5 * M_PI_SQ - 4.0;
    if (x == 2.0)
        return M_PI_SQ / 6.0 - 1.0;

    /// Use integral representation
    return -integral([x] (double t)
    {
        if (t >= 1)
            return -1.0;
        if (t <= 0)
            return 0.0;
        double logT = std::log(t);
        double y = x * logT;
        y = std::exp(y);
        y *= logT;
        y /= (1.0 - t);
        return y;
    }, 0, 1) + 1.0 / (x * x);
}

long double lowerIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return 0.0;
    if (a == 1)
        return 1.0 - std::exp(-x);
    if (a == 0.5)
        return M_SQRTPI * std::erf(std::sqrt(x));
    if (x > 1.1 && x > a) {
        return std::tgamma(a) - upperIncGamma(a, x);
    }
    return std::exp(logLowerIncGamma(a, x));
}

long double logLowerIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return -INFINITY;
    if (a == 1)
        return std::log1p(-std::exp(-x));
    if (a == 0.5)
        return 0.5 * M_LNPI + std::log(std::erf(std::sqrt(x)));

    if (x > 1.1 && x > a) {
        double y = std::tgamma(a) - upperIncGamma(a, x);
        return std::log(y);
    }

    long double sum = 0;
    long double term = 1.0 / a;
    double n = a + 1;
    while (std::fabs(term) > MIN_POSITIVE)
    {
        sum = sum + term;
        term *= x / n;
        ++n;
    }
    return a * std::log(x) - x + std::log(sum);
}

long double upperIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return std::tgamma(x);
    if (a == 0.5)
        return M_SQRTPI * std::erfc(std::sqrt(x));
    if (a == 1)
        return std::exp(-x);
    if (x <= 1.1 || x <= a) {
        return std::tgamma(a) - lowerIncGamma(a, x);
    }
    return std::exp(logUpperIncGamma(a, x));
}

long double logUpperIncGamma(double a, double x)
{
    if (x < 0)
        return NAN;
    if (x == 0)
        return std::lgamma(x);
    if (a == 1)
        return -x;
    if (a == 0.5)
        return 0.5 * M_LNPI + std::log(std::erfc(std::sqrt(x)));
    if (std::max(a, x) < 1e-5)
        return a * std::log(x);
    if (x > 1e4) {
        return (a - 1) * std::log(x) - x;
    }

    if (x <= 1.1 || x <= a) {
        double y = std::tgamma(a) - lowerIncGamma(a, x);
        return std::log(y);
    }

    /// the modified Lentz's method
    double fraction = 1.0 + x - a;
    double C = fraction, D = 0, mult = 1;
    double xmap1 = x - a + 1.0;
    int i = 1;
    do {
        double s = i * (a - i);
        double b = (i << 1) + xmap1;
        D = b + s * D;
        C = b + s / C;
        D = 1.0 / D;
        mult = C * D;
        fraction *= mult;
    } while (std::fabs(mult - 1.0) > MIN_POSITIVE && ++i < 10000);
    return a * std::log(x) - x - std::log(fraction);
}

double betaFun(double a, double b)
{
    if (a <= 0 || b <= 0)
        return NAN;
    double lgammaA = std::lgamma(a);
    double lgammaB = (a == b) ? lgammaA : std::lgamma(b);
    return std::exp(lgammaA + lgammaB - std::lgamma(a + b));
}

double regularizedBetaFun(double x, double a, double b)
{
    if (a <= 0 || b < 0 || x < 0.0 || x > 1.0)
        return NAN;
    if (b == 0.0 || x == 0.0)
        return 0.0;
    if (x == 1.0)
        return 1.0;
    return incompleteBetaFun(x, a, b) / betaFun(a, b);
}

double incompleteBetaFun(double x, double a, double b)
{
    if (a <= 0 || b < 0 || x < 0.0 || x > 1.0) /// if incorrect parameters
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return (b == 0) ? NAN : betaFun(a, b);
    if (a < 1)
    {
        double y = incompleteBetaFun(x, a + 1, b) * (a + b);
        double z = a * std::log(x) + b * std::log1p(-x);
        y += std::exp(z);
        return y / a;
    }
    if (b < 1)
    {
        if (b == 0) /// series expansion
        {
            double denom = a, numen = 1;
            double sum = 1.0 / a, add = 1;
            do {
                numen *= x;
                add = numen / (++denom);
                sum += add;
            } while (add > MIN_POSITIVE);
            return std::pow(x, a) * sum;
        }
        double y = incompleteBetaFun(x, a, b + 1) * (a + b);
        double z = a * std::log(x) + b * std::log1p(-x);
        y -= std::exp(z);
        return y / b;
    }

    double minBound = 0, maxBound = x;
    bool invert = false;
    /// if x > mode
    if (x > (a - 1) / (a + b - 2)) {
        maxBound = 1;
        minBound = x;
        invert = true;
    }
    double y = 0;

    if (a != b)
    {
        y = integral([a, b] (double t)
        {
            double z = (a - 1) * std::log(t) + (b - 1) * std::log1p(-t);
            return std::exp(z);
        },
        minBound, maxBound);
    }
    else {
        y = integral([a, b] (double t)
        {
            return std::pow(t - t * t, a - 1);
        },
        minBound, maxBound);
    }
    return (invert) ? betaFun(a, b) - y : y;
}

long double gammaHalf(size_t k)
{
    if (k & 1)
    {
        size_t n = (k - 1) >> 1;
        long double res = factorial(k - 1);
        res /= (factorial(n) * (1 << (n << 1)));
        return res * M_SQRTPI;
    }

    return factorial((k >> 1) - 1);
}

/**
 * @brief adaptiveSimpsonsAux
 * auxiliary function for calculation of integral
 * @param funPtr
 * @param a lower boundary
 * @param b upper boundary
 * @param epsilon
 * @param S
 * @param fa
 * @param fb
 * @param fc
 * @param bottom
 * @return
 */
long double adaptiveSimpsonsAux(const std::function<double (double)> &funPtr, double a, double b,
                                          double epsilon, double S, double fa, double fb, double fc, int bottom)
{
    double c = .5 * (a + b), h = (b - a) / 12.0;
    double d = .5 * (a + c), e = .5 * (c + b);
    double fd = funPtr(d), fe = funPtr(e);
    double Sleft = h * (fa + 4 * fd + fc);
    double Sright = h * (fc + 4 * fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || std::fabs(S2 - S) <= 15.0 * epsilon)
        return S2 + (S2 - S) / 15.0;
    epsilon *= .5;
    --bottom;

    return adaptiveSimpsonsAux(funPtr, a, c, epsilon, Sleft, fa, fc, fd, bottom) +
           adaptiveSimpsonsAux(funPtr, c, b, epsilon, Sright, fc, fb, fe, bottom);
}

long double integral(const std::function<double (double)> &funPtr,
                               double a, double b, double epsilon, int maxRecursionDepth)
{
    if (a > b)
        std::swap(a, b);
    if (a == b)
        return 0.0;
    double c = .5 * (a + b), h = (b - a) / 6.0;
    double fa = funPtr(a), fb = funPtr(b), fc = funPtr(c);
    double S = h * (fa + 4 * fc + fb);
    return adaptiveSimpsonsAux(funPtr, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

bool findRoot(const std::function<DoubleTriplet (double)> &funPtr, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);
    static constexpr int maxIter = 1e5;
    int iter = 0;
    double step = epsilon + 1;
    DoubleTriplet y = funPtr(root);
    double f, fx, fxx;
    std::tie(f, fx, fxx) = y;
    if (std::fabs(f) < epsilon)
        return root;
    do {
        double alpha = 1.0;
        double oldRoot = root;
        double oldFun = f;
        double numerator = 2 * f * fx;
        double denominator = 2 * fx * fx - f * fxx;
        step = numerator / denominator;
        do {
            root = oldRoot - alpha * step;
            y = funPtr(root);
            // TODO: make bounds for fx and fxx
            std::tie(f, fx, fxx) = y;
            if (std::fabs(f) < epsilon)
                return true;
            alpha *= 0.5;
        } while ((std::fabs(fx) <= epsilon || std::fabs(oldFun) < std::fabs(f)) && alpha > 0);
    } while (std::fabs(step) > epsilon && ++iter < maxIter);

    return (iter == maxIter) ? false : true;
}

bool findRoot(const std::function<DoublePair (double)> &funPtr, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);
    static constexpr int MAX_ITER = 1e5;
    static constexpr double MAX_GRAD = 1e4;
    int iter = 0;
    double step = epsilon + 1;
    DoublePair y = funPtr(root);
    double fun = y.first;
    double grad = std::min(MAX_GRAD, std::max(-MAX_GRAD, y.second));
    if (std::fabs(fun) < epsilon)
        return root;
    do {
        double alpha = 1.0;
        double oldRoot = root;
        double oldFun = fun;
        step = fun / grad;
        do {
            root = oldRoot - alpha * step;
            y = funPtr(root);
            fun = y.first;
            grad = std::min(MAX_GRAD, std::max(-MAX_GRAD, y.second));
            if (std::fabs(fun) < epsilon)
                return true;
            alpha *= 0.5;
        } while ((std::fabs(grad) <= epsilon || std::fabs(oldFun) < std::fabs(fun)) && alpha > 0);
    } while (std::fabs(step) > epsilon && ++iter < MAX_ITER);

    return (iter == MAX_ITER) ? false : true;
}

bool findRoot(const std::function<double (double)> &funPtr, double a, double b, double &root, double epsilon)
{
    /// Sanity check
    epsilon = std::max(epsilon, MIN_POSITIVE);

    double fa = funPtr(a);
    if (fa == 0) {
        root = a;
        return true;
    }

    double fb = funPtr(b);
    if (fb == 0) {
        root = b;
        return true;
    }

    if (fa * fb > 0) {
        /// error - the root is not bracketed
        return false;
    }
    if (std::fabs(fa) < std::fabs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    double c = a, fc = fa;
    bool mflag = true;
    double s = b, fs = 1, d = 0;
    while (std::fabs(b - a) > epsilon) {
        if (!areClose(fc, fa) && !areClose(fb, fc))
        {
            /// inverse quadratic interpolation
            double numerator = a * fb * fc;
            double denominator = (fa - fb) * (fa - fc);
            s = numerator / denominator;

            numerator = b * fa * fc;
            denominator = (fb - fa) * (fb - fc);
            s += numerator / denominator;

            numerator = c * fa * fb;
            denominator = (fc - fa) * (fc - fb);
            s += numerator / denominator;
        }
        else
        {
            /// secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        double absDiffSB2 = std::fabs(s - b);
        absDiffSB2 += absDiffSB2;
        double absDiffBC = std::fabs(b - c);
        double absDiffCD = std::fabs(c - d);
        if (s < 0.25 * (3 * a + b) || s > b ||
            (mflag && absDiffSB2 >= absDiffBC) ||
            (!mflag && absDiffSB2 >= absDiffCD) ||
            (mflag && absDiffBC < epsilon) ||
            (!mflag && absDiffCD < epsilon))
        {
            s = 0.5 * (a + b);
            mflag = true;
        }
        else {
            mflag = false;
        }

        fs = funPtr(s);
        if (std::fabs(fs) < epsilon) {
            root = s;
            return true;
        }

        d = c;
        c = b;
        fc = fb;

        if (fa * fs < 0) {
            b = s;
            fb = fs;
        }
        else {
            a = s;
            fa = fs;
        }

        if (std::fabs(fa) < std::fabs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    root = (std::fabs(fs) < std::fabs(fb)) ? s : b;
    return true;
}

/**
 * @brief parabolicMinimum
 * @param a < b < c
 * @param fa f(a)
 * @param fb f(b)
 * @param fc f(c)
 * @return minimum of interpolated parabola
 */
double parabolicMinimum(double a, double b, double c, double fa, double fb, double fc)
{
    double bma = b - a, cmb = c - b;
    double aux1 = bma * (fb - fc);
    double aux2 = cmb * (fb - fa);
    double numerator = bma * aux1 - cmb * aux2;
    double denominator = aux1 + aux2;
    return b - 0.5 * numerator / denominator;
}

/**
 * @brief findBounds
 * Search of segment containing minimum of function
 * @param funPtr mapping x |-> f(x)
 * @param abc such points, that a < b < c, f(a) > f(b) and f(c) > f(b)
 * @param fabc values of a, b and c
 * @param startPoint
 * @return true when segment is found
 */
bool findBounds(const std::function<double (double)> &funPtr, DoubleTriplet &abc, DoubleTriplet &fabc, double startPoint)
{
    static constexpr double K = 0.5 * (M_SQRT5 + 1);
    static constexpr int L = 100;
    double a = startPoint, fa = funPtr(a);
    double b = a + 1.0, fb = funPtr(b);
    double c, fc;
    if (fb < fa) {
        c = b + K * (b - a);
        fc = funPtr(c);
        /// we go to the right
        while (fc < fb) {
            /// parabolic interpolation
            double u = parabolicMinimum(a, b, c, fa, fb, fc);
            double cmb = c - b;
            double fu, uLim = c + L * cmb;
            if (u < c && u > b) {
                fu = funPtr(u);
                if (fu < fc) {
                    abc = std::make_tuple(b, u, c);
                    fabc = std::make_tuple(fb, fu, fc);
                    return true;
                }
                if (fu > fb) {
                    abc = std::make_tuple(a, b, u);
                    fabc = std::make_tuple(fa, fb, fu);
                    return true;
                }
                u = c + K * cmb;
                fu = funPtr(u);
            }
            else if (u > c && u < uLim) {
                fu = funPtr(u);
                if (fu < fc) {
                    b = c; c = u; u = c + K * cmb;
                    fb = fc, fc = fu, fu = funPtr(u);
                }
            }
            else if (u > uLim) {
                u = uLim;
                fu = funPtr(u);
            }
            else {
                u = c + K * cmb;
                fu = funPtr(u);
            }
            a = b; b = c; c = u;
            fa = fb; fb = fc; fc = fu;
        }
        abc = std::make_tuple(a, b, c);
        fabc = std::make_tuple(fa, fb, fc);
        return true;
    }
    else {
        c = b; fc = fb;
        b = a; fb = fa;
        a = b - K * (c - b);
        fa = funPtr(a);
        /// go to the left
        while (fa < fb) {
            /// parabolic interpolation
            double u = parabolicMinimum(a, b, c, fa, fb, fc);
            double bma = b - a;
            double fu, uLim = a - L * bma;
            if (u < b && u > a) {
                fu = funPtr(u);
                if (fu < fa) {
                    abc = std::make_tuple(a, u, b);
                    fabc = std::make_tuple(fa, fu, fb);
                    return true;
                }
                if (fu > fb) {
                    abc = std::make_tuple(u, b, c);
                    fabc = std::make_tuple(fu, fb, fc);
                    return true;
                }
                u = a - K * bma;
                fu = funPtr(u);
            }
            else if (u < a && u > uLim) {
                fu = funPtr(u);
                if (fu < fa) {
                    b = a; a = u; u = a - K * bma;
                    fb = fa, fa = fu, fu = funPtr(u);
                }
            }
            else if (u < uLim) {
                u = uLim;
                fu = funPtr(u);
            }
            else {
                u = a - K * bma;
                fu = funPtr(u);
            }
            c = b; b = a; a = u;
            fc = fb; fb = fa; fa = fu;
        }
        abc = std::make_tuple(a, b, c);
        fabc = std::make_tuple(fa, fb, fc);
        return true;
    }
}

bool findMin(const std::function<double (double)> &funPtr, DoubleTriplet abc, DoubleTriplet fabc, double &root, double epsilon)
{
    static constexpr double K = 0.5 * (3 - M_SQRT5);
    double a, x, c;
    std::tie(a, x, c) = abc;
    double fa, fx, fc;
    std::tie(fa, fx, fc) = fabc;
    double w = x, v = x, fw = fx, fv = fx;
    double d = c - a, e = d;
    double u = a - 1;
    do {
        double g = e;
        e = d;
        bool acceptParabolicU = false;
        if (x != w && x != v && w != v &&
            fx != fw && fx != fv && fw != fv) {
            if (v < w) {
                if (x < v)
                    u = parabolicMinimum(x, v, w, fx, fv, fw);
                else if (x < w)
                    u = parabolicMinimum(v, x, w, fv, fx, fw);
                else
                    u = parabolicMinimum(v, w, x, fv, fw, fx);
            }
            else {
                if (x < w)
                    u = parabolicMinimum(x, w, v, fx, fv, fw);
                else if (x < v)
                    u = parabolicMinimum(w, x, v, fw, fx, fv);
                else
                    u = parabolicMinimum(w, v, x, fw, fv, fx);
            }
            double absumx = std::fabs(u - x);
            if (u >= a + epsilon && u <= c - epsilon && absumx < 0.5 * g) {
                acceptParabolicU = true; /// accept u
                d = absumx;
            }
        }

        if (!acceptParabolicU) {
            /// use golden ratio instead of parabolic approximation
            if (x < 0.5 * (c + a)) {
                d = c - x;
                u = x + K * d; /// golden ratio [x, c]
            }
            else {
                d = x - a;
                u = x - K * d; /// golden ratio [a, x]
            }
        }

        if (std::fabs(u - x) < epsilon) {
            u = x + epsilon * sign(u - x); /// setting the closest distance between u and x
        }

        double fu = funPtr(u);
        if (fu <= fx) {
            if (u >= x)
                a = x;
            else
                c = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        }
        else {
            if (u >= x)
                c = u;
            else
                a = u;
            if (fu <= fw || w == x) {
                v = w; w = u;
                fv = fw; fw = fu;
            }
            else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    } while (0.49 * (c - a) > epsilon);
    root = x;
    return true;
}

bool findMin(const std::function<double (double)> &funPtr, double closePoint, double &root, double epsilon)
{
    DoubleTriplet abc, fabc;
    if (!findBounds(funPtr, abc, fabc, closePoint))
        return false;
    return findMin(funPtr, abc, fabc, root, epsilon);
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
    double res = 1.0;
    for (int i = 2; i <= number; ++i)
        res += std::pow(i, -exponent);
    return res;
}

double modifiedBesselFirstKind(double x, double n)
{
    double roundN = std::round(n);
    bool nIsInt = areClose(n, roundN);

    if (x < 0) {
        if (nIsInt)
        {
            int nInt = roundN;
            return (nInt % 2) ? -modifiedBesselFirstKind(-x, n) : modifiedBesselFirstKind(-x, n);
        }
        else
            return 0.0;
    }

    if (x == 0) {
        if (n == 0)
            return 1.0;
        return (n > 0 || nIsInt) ? 0.0 : INFINITY;
    }

    if (n == 0.5)
        return std::sqrt(M_2_PI / x) * std::sinh(x);

    if (n == -0.5)
        return std::sqrt(M_2_PI / x) * std::cosh(x);

    /// small x
    if (x < 10)
    {
        double halfX = 0.5 * x;
        double halfXSq = halfX * halfX;
        double addon = 1.0;
        double sum = 1.0;
        double i = 1.0;
        while (std::fabs(addon) > MIN_POSITIVE) {
            addon *= halfXSq;
            addon /= i * (i + n);
            ++i;
            sum += addon;
        }
        double y = std::log(sum);
        y += n * std::log(halfX);
        y -= std::lgamma(n + 1);
        return std::exp(y);
    }

    // large x - divergent sequence!!
    double addon = 1.0;
    double sum = addon;
    double denominator = 0.125 / x;
    double n2Sq = 4.0 * n * n;
    double i = 1.0, j = 1.0;
    while (std::fabs(addon) > MIN_POSITIVE) {
        double numerator = j * j - n2Sq;
        double frac = numerator * denominator / i;
        if (frac > 1)
            break;
        addon *= frac;
        sum += addon;
        ++i;
        j += 2;
    }
    double y = M_PI * x;
    y = std::sqrt(y + y);
    y = std::exp(x) / y;
    return y * sum;
}

double modifiedBesselSecondKind(double x, double n)
{
    if (n == 0.5)
        return M_SQRTPI * M_SQRT1_2 * std::exp(-x) / std::sqrt(x);

    double y = modifiedBesselFirstKind(x, -n);
    y -= modifiedBesselFirstKind(x, n);
    y /= std::sin(n * M_PI);
    return 0.5 * M_PI * y;
}

double zetaRiemann(double s)
{
    if (s == 1)
        return INFINITY;
    // TODO!! http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf
    if (s < 1)
        return NAN;
    int N = 100;
    double y = harmonicNumber(s, N);
    double NS = std::pow(N, -s);
    return y + N * NS / (s - 1) + 0.5 * NS;
}

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

}

