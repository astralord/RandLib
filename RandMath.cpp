#include "RandMath.h"

RandMath::RandMath()
{
}

long double RandMath::stirlingFactorial(int n)
{
    double fact = M_SQRT2PI * std::sqrt(static_cast<double>(n));
    return fact * std::pow(n / M_E, n);
}

long double RandMath::getTenthFactorial(int n)
{
    switch (n) {
    case 0:
        return 1;
    case 10:
        return 3628800;
    case 20:
        return 2432902008176640000;
    case 30:
        return 265252859812191058636308480000000.;
    case 40:
        return 815915283247897734345611269596115894272000000000.;
    case 50:
        return 30414093201713378043612608166064768844377641568960512000000000000.;
    case 60:
        return 8320987112741390144276341183223364380754172606361245952449277696409600000000000000.;
    case 70:
        return 11978571669969891796072783721689098736458938142546425857555362864628009582789845319680000000000000000.;
    case 80:
        return 71569457046263802294811533723186532165584657342365752577109445058227039255480148842668944867280814080000000000000000000.;
    case 90:
        return 1485715964481761497309522733620825737885569961284688766942216863704985393094065876545992131370884059645617234469978112000000000000000000000.;
    case 100:
        return 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000.;
    case 110:
        return 15882455415227429404253703127090772871724410234473563207581748318444567162948183030959960131517678520479243672638179990208521148623422266876757623911219200000000000000000000000000.;
    case 120:
        return 6689502913449127057588118054090372586752746333138029810295671352301633557244962989366874165271984981308157637893214090552534408589408121859898481114389650005964960521256960000000000000000000000000000.;
    case 130:
        return 6466855489220473672507304395536485253155359447828049608975952322944781961185526165512707047229268452925683969240398027149120740074042105844737747799459310029635780991774612983803150965145600000000000000000000000000000000.;
    case 140:
        return 13462012475717524605876073858941615558355851148193967190051391468057460367090535696797920946629681836680869097041958983702264048370902871114013579941370766400374327741701139895604871545254810788060989321379840000000000000000000000000000000000.;
    case 150:
        return 57133839564458545904789328652610540031895535786011264182548375833179829124845398393126574488675311145377107878746854204162666250198684504466355949195922066574942592095735778929325357290444962472405416790722118445437122269675520000000000000000000000000000000000000.;
    case 160:
        return 471472363599206132240694321176194377951192623045460204976904578317542573467421580346978030238114995699562728104819596262106947389303901748942909887857509625114880781313585012959529941660203611234871833992565791817698209861793313332044734813700096000000000000000000000000000000000000000.;
    case 170:
        return 7257415615307998967396728211129263114716991681296451376543577798900561843401706157852350749242617459511490991237838520776666022565442753025328900773207510902400430280058295603966612599658257104398558294257568966313439612262571094946806711205568880457193340212661452800000000000000000000000000000000000000000.;
    case 180:
        return 200896062499134299656951336898466838917540340798867777940435335160044860953395980941180138112097309735631594101037399609671032132186331495273609598531966730972945653558819806475064353856858157445040809209560358463319644664891114256430017824141796753818192338642302693327818731986039603200000000000000000000000000000000000000000000.l;
    case 190:
        return 9680322675255249156123346514615331205418161260462873360750859919944104623425228207640470674933540169424682360525991982916161596983449594045525553704253602287443197783274656957056546338783001340434094795097553229620273057440272298773179365935914105128629426348958748638226084106818484328004851174161755668480000000000000000000000000000000000000000000000.l;
    case 200:
        return 788657867364790503552363213932185062295135977687173263294742533244359449963403342920304284011984623904177212138919638830257642790242637105061926624952829931113462857270763317237396988943922445621451664240254033291864131227428294853277524242407573903240321257405579568660226031904170324062351700858796178922222789623703897374720000000000000000000000000000000000000000000000000.l;
    case 210:
        return 105823620292236563784274284243348353057589905787169019562352737522144487532400210147849369011714673954768265316577892528273760626189481169051055226066650741189573897273684791411180134039439160066561895838501000817711682625725670477616267598661259194975646029749546282594356217374097544153589482020891750774735012558313460846824864172030239122128896000000000000000000000000000000000000000000000000000.l;
    case 220:
        return 22838603359146414573972658651153337270429730715462287017736347161260276926030248458777765497919211029457065581960747795750095505232241970499561769723020565876672261660609763234049775547325430135571331468257475537994508495233770658945310210552725163342784668756149049213658078338458534285571551800849578848226429898670032945513859929938621783523490272646966918544936140800000000000000000000000000000000000000000000000000000.l;
    case 230:
        return 7758587304686725201813174298892781442413952130995533365303964524344944412641389739603152000644515957408814002319492032321234250506968028455594445689972313374305301019340949789291189972149450405025159624155827152329676580440959428615802893638146558163235483142136540783687811997927615346859658417205832954125915861983307177232587595821512723429698627780530255874167602077755356592824804966400000000000000000000000000000000000000000000000000000000.l;
    case 240:
        return 4067885363647058120493575921486885310172051259182827146069755969081486918925585104009100729728348522923820890245870098659147156051905732563147381599098459244752463027688115705371704628286326621238456543307267608612545168337779669138759451760395968217423617954330737034164596496963986817722252221059768080852489940995605579171999666916004042965293896799800598079985264195119506681577622056215044851618236292136960000000000000000000000000000000000000000000000000000000000.l;
    case 250:
        return 3232856260909107732320814552024368470994843717673780666747942427112823747555111209488817915371028199450928507353189432926730931712808990822791030279071281921676527240189264733218041186261006832925365133678939089569935713530175040513178760077247933065402339006164825552248819436572586057399222641254832982204849137721776650641276858807153128978777672951913990844377478702589172973255150283241787320658188482062478582659808848825548800000000000000000000000000000000000000000000000000000000000000.l;
    default:
        return 0;
    }
}

long double RandMath::fastFactorial(int n)
{
    if (n < 0)
        return 0;

    if (n > 255)
        return stirlingFactorial(n);

    int residue = n % 10;
    if (n <= 5)
    {
        /// go up
        int nPrev = n - residue;
        double fact = getTenthFactorial(nPrev);
        for (int i = 1; i <= residue; ++i)
            fact *= nPrev + i;
        return fact;
    }
    else
    {
        /// go  down
        int nNext = n - residue + 10;
        double denom = 1;
        for (int i = 0; i < 10 - residue; ++i)
            denom *= nNext - i;
        return getTenthFactorial(nNext) / denom;
    }
}

long double RandMath::doubleFactorial(int n)
{
    long double n_fact = fastFactorial(n);
    if (n % 2 == 0)
        return (1 << n) * n_fact;
    return fastFactorial(2 * n + 1) / (2 * n * n_fact);
}

long double RandMath::binomialCoef(int n, int k)
{
    long double n_fact = fastFactorial(n);
    long double k_fact = fastFactorial(k);
    long double k_n_fact = fastFactorial(n - k);
    return n_fact / (k_fact * k_n_fact);
}

long double RandMath::lowerIncGamma(double a, double x)
{
    double sum = 0;
    double term = 1.0 / a;
    int n = 1;
    while (std::fabs(term) > MIN_POSITIVE)
    {
        sum = sum + term;
        term *= (x / (a + n));
        ++n;
    }
    return std::pow(x, a) * std::exp(-x) * sum;
}

long double RandMath::upperIncGamma(double a, double x)
{
    // TODO: find useful approximation
    return std::tgamma(a) - lowerIncGamma(a, x);
}

long double RandMath::betaFun(double x, double y)
{
    if (x > y)
    {
        long double res = std::tgamma(x);
        res /= std::tgamma(x + y);
        return res * std::tgamma(y);
    }
    else
    {
        long double res = std::tgamma(y);
        res /= std::tgamma(x + y);
        return res * std::tgamma(x);
    }
}

long double RandMath::gammaHalf(int k)
{
    if (k < 0)
        return 0;

    if (k % 2 == 0)
        return fastFactorial((k >> 1) - 1);

    int n = (k - 1) >> 1;
    long double res = fastFactorial(k - 1);
    res /= (fastFactorial(n) * (1 << (n + n)));
    return res * M_SQRTPI;
}

long double RandMath::erfInv(double p)
{
    if (p < 0.5)
        return -erfInv(1 - p);
    if (p <= 0)
        return -INFINITY;
    if (p >= 1)
        return INFINITY;
    double t = M_SQRT2 * std::sqrt(-std::log(p));
    static constexpr double c[] = {2.515517, 0.802853, 0.010328};
    static constexpr double d[] = {1.432788, 0.189269, 0.001308};
    long double numen = (c[2] * t + c[1]) * t + c[0];
    long double denom = ((d[2] * t + d[1]) * t + d[0]) * t + 1.0;
    return t - numen / denom;
}

long double RandMath::erfcinv(double p)
{
    return erfInv(1 - p);
}

/*
long double RandMath::adaptiveSimpsonsAux(std::function<double (const RandomVariable &, double)> fun, const RandomVariable &rv,
                               double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom)
{
    // TODO: rewrite recursion into loop
    double c = .5 * (a + b), h = (b - a) / 12.0;
    double d = .5 * (a + c), e = .5 * (c + b);
    double fd = fun(rv, d), fe = fun(rv, e);
    double Sleft = h * (fa + 4 * fd + fc);
    double Sright = h * (fc + 4 * fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || std::fabs(S2 - S) <= 15.0 * epsilon)
        return S2 + (S2 - S) / 15.0;
    epsilon *= .5;
    --bottom;

    return adaptiveSimpsonsAux(fun, rv, a, c, epsilon, Sleft, fa, fc, fd, bottom) +
    adaptiveSimpsonsAux(fun, rv, c, b, epsilon, Sright, fc, fb, fe, bottom);
}

long double RandMath::integral(std::function<double (const RandomVariable &, double)> fun, const RandomVariable &rv,
                               double a, double b, double epsilon, int maxRecursionDepth)
{
    double c = .5 * (a + b), h = (b - a) / 6.0;
    double fa = fun(rv, a), fb = fun(rv, b), fc = fun(rv, c);
    double S = h * (fa + 4 * fc + fb);
    return adaptiveSimpsonsAux(fun, rv, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}
*/