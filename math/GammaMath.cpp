#include "GammaMath.h"

namespace RandMath
{

/**
 * @fn FACTORIAL_TABLESIZE maximum value for input parameter to use table method
 */
constexpr int FACTORIAL_TABLESIZE = 255;

/**
 * @fn FACTORIAL_TABLE (n * 10)! for n from 0 to 25
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
 * @fn factorialForSmallValue
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
    long double nFact = factorial(n);
    if (n & 1) {
        n <<= 1;
        return factorial(n + 1) / (n * nFact);
    }
    return (1 << n) * nFact;
}

constexpr long double LOGFACTORIAL_TABLE[] =
{
    0.l,
    0.l,
    0.693147180559945286227l,
    1.79175946922805495731l,
    3.17805383034794575181l,
    4.78749174278204581157l,
    6.57925121201010121297l,
    8.52516136106541466688l,
    10.6046029027452508586l,
    12.8018274800814690906l,
    15.1044125730755158799l,
    17.5023078458738865493l,
    19.9872144956618846834l,
    22.5521638531234209779l,
    25.1912211827386798291l,
    27.8992713838408903371l,
    30.6718601060806754788l,
    33.5050734501368836504l,
    36.3954452080330526087l,
    39.3398841871994946473l,
    42.3356164607534850575l,
    45.3801388984769147328l,
    48.4711813518352201413l,
    51.6066755677643769218l,
    54.7847293981123115714l,
    58.0036052229805250136l,
    61.2617017610020013763l,
    64.5575386270063233951l,
    67.889743137181540078l,
    71.2570389671680146648l,
    74.6582363488301581356l,
    78.0922235533153070719l,
    81.5579594561150429399l,
    85.0544670175815298307l,
    88.5808275421976674124l,
    92.1361756036870929165l,
    95.719694542143201943l,
    99.3306124547874276232l,
    102.96819861451380973l,
    106.63176026064346047l,
    110.320639714757405159l,
    114.034211781461706892l,
    117.77188139974506953l,
    121.533081515438624365l,
    125.317271149356884052l,
    129.123933639127216111l,
    132.952575035616291643l,
    136.802722637326382937l,
    140.67392364823425055l,
    144.565743946344895221l,
    148.47776695177302031l,
    152.409592584497374901l,
    156.360836303078798437l,
    160.331128216630901306l,
    164.320112263195198921l,
    168.32744544842768164l,
    172.352797139162817075l,
    176.395848406997345137l,
    180.456291417543781108l,
    184.533828861449478609l,
    188.628173423671597675l,
    192.739047287844897483l,
    196.866181672889979382l,
    201.009316399281516397l,
    205.168199482641171016l,
    209.342586752536817585l,
    213.532241494563237438l,
    217.73693411395424846l,
    221.956441819130361637l,
    226.190548323727625757l,
    230.439043565776955802l,
    234.701723442818263266l,
    238.978389561834319466l,
    243.268849002982733509l,
    247.57291409618684952l,
    251.890402209723191618l,
    256.221135550009535109l,
    260.564940971863222785l,
    264.921649798552778066l,
    269.29109765101981111l,
    273.673124285693745605l,
    278.067573440366174964l,
    282.474292687630452292l,
    286.893133295426991936l,
    291.32395009427034438l,
    295.766601350760595324l,
    300.220948647014154176l,
    304.686856765668721891l,
    309.164193580146900331l,
    313.652829949879048854l,
    318.152639620209299665l,
    322.663499126726208033l,
    327.18528770377525916l,
    331.717887196928472804l,
    336.261181979198454428l,
    340.815058870799020951l,
    345.379407062266864159l,
    349.954118040770254083l,
    354.539085519440789085l,
    359.13420536957545437l,
    363.739375555563526632l,
    368.354496072404685947l,
    372.979468885689016133l,
    377.614197873918612913l,
    382.258588773060012045l,
    386.912549123217502256l,
    391.575988217329552299l,
    396.248817051791547783l,
    400.930948278915707306l,
    405.622296161144902271l,
    410.322776526937275321l,
    415.032306728249579919l,
    419.750805599544776214l,
    424.478193418257092162l,
    429.214391866651567398l,
    433.959323995014813136l,
    438.712914186121224702l,
    443.475088120918996992l,
    448.245772745384613245l,
    453.024896238496125989l,
    457.81238798127822065l,
    462.608178526874951331l,
    467.412199571608141468l,
    472.224383926980635806l,
    477.0446654925856933l,
    481.87297922988790333l,
    486.70926113683941594l,
    491.553448223297948516l,
    496.405478487217578731l,
    501.265290891579240906l,
    506.132825342034834648l,
    511.008022665236012472l,
    515.890824587822407921l,
    520.781173716044122557l,
    525.67901351599505233l,
    530.584288294433576993l,
    535.496943180169523657l,
    540.41692410599773666l,
    545.344177791154834267l,
    550.278651724285509772l,
    555.220294146894843834l,
    560.169054037273099311l,
    565.124881094874240262l,
    570.087725725134191634l,
    575.057539024710195008l,
    580.034272767130801185l,
    585.017879388839105559l,
    590.008311975617743883l,
    595.005524249382006019l,
    600.009470555327425245l,
    605.020105849423657673l,
    610.037385686238621929l,
    615.061266207084827329l,
    620.091704128477431368l,
    625.128656730890952531l,
    630.172081847810204636l,
    635.221937855059650246l,
    640.278183660408103606l,
    645.34077869343502698l,
    650.409682895655237189l,
    655.484856710889062015l,
    660.566261075873512709l,
    665.653857411105946085l,
    670.747607611912712855l,
    675.847474039736766827l,
    680.953419513637527416l,
    686.065407301994014233l,
    691.183401114410685295l,
    696.307365093814041757l,
    701.437263808737156978l,
    706.573062245787468783l,
    711.714725802289990497l,
    716.862220279103553366l,
    722.015511873601212756l,
    727.174567172815841332l,
    732.339353146739199474l,
    737.509837141777325087l,
    742.685986874351215192l,
    747.867770424643367733l,
    753.055156230484158186l,
    758.248113081374413014l,
    763.446610112640087209l,
    768.650616799716999594l,
    773.860102952558349898l,
    779.075038710167291356l,
    784.295394535245691259l,
    789.521141208958852076l,
    794.752249825813464668l,
    799.988691788643450309l,
    805.230438803703009398l,
    810.477462875863579939l,
    815.729736303910158313l,
    820.987231675938005537l,
    826.249921864842804098l,
    831.517780023906198039l,
    836.790779582469781417l,
    842.068894241700377279l,
    847.352097970438421726l,
    852.640365001132863654l,
    857.933669825857350588l,
    863.231987192405426867l,
    868.535292100464630494l,
    873.843559797865736982l,
    879.156765776907491272l,
    884.474885770751825476l,
    889.797895749890244588l,
    895.125771918679674855l,
    900.458490711945046314l,
    905.796028791646449463l,
    911.138363043611207104l,
    916.485470574328701332l,
    921.837328707804772421l,
    927.193914982476826481l,
    932.555207148186127597l,
    937.921183163207956568l,
    943.291821191335770891l,
    948.667099599019934431l,
    954.046996952560334648l,
    959.431492015349476787l,
    964.820563745166055014l,
    970.214191291518204707l,
    975.612353993035981148l,
    981.01503137490840345l,
    986.422203146368360649l,
    991.833849198223447274l,
    997.249949600427953555l,
    1002.67048459970021668l,
    1008.0954346171816951l,
    1013.52478024613617436l,
    1018.9585022496903548l,
    1024.39658155861343403l,
    1029.83899926913522904l,
    1035.28573664080158778l,
    1040.73677509436743094l,
    1046.19209620972492303l,
    1051.65168172386916012l,
    1057.11551352889478039l,
    1062.58357367002986393l,
    1068.05584434370143754l,
    1073.53230789563281178l,
    1079.01294681897479677l,
    1084.49774375246556701l,
    1089.98668147862213118l,
    1095.47974292196272472l,
    1100.97691114725603256l,
    1106.47816935780065251l,
    1111.9835008937332077l,
    1117.49288923036101551l,
    1123.00631797652590649l,
    1128.52377087299055347l,
    1134.04523179085276752l,
    1139.57068472998480502l,
    1145.10011381749609427l,
    1150.63350330622370166l,
    1156.17083757324212456l,
    1161.71210111840059653l
};

constexpr size_t LOGFACTORIAL_TABLESIZE = 255;

long double lfact(size_t n)
{
    return (n > LOGFACTORIAL_TABLESIZE) ? std::lgammal(n + 1) : LOGFACTORIAL_TABLE[n];
}

long double ldfact(size_t n)
{
    if (n & 1) {
        return lfact(n) - ldfact(n - 1);
    }
    size_t k = n >> 1;
    return k * M_LN2 + lfact(k);
}

long double binom(size_t n, size_t k)
{
    /// first check trivial cases
    if (k == 0 || k == n)
        return 1.0;
    if (k == 1 || k == n - 1)
        return n;
    if (k == 2 || k == n - 2)
        return 0.5 * n * (n - 1);
    /// next run general procedure
    double lfactn = lfact(n);
    double lfactk = lfact(k);
    double lfactnmk = lfact(n - k);
    return std::exp(lfactn - lfactk - lfactnmk);
}

/**
 * @fn digammamLogForLargeX
 * @param x
 * @return digamma(x) - log(x) for x > 7
 */
double digammamLogForLargeX(double x)
{
    /// set up tables
    static constexpr long double taylorCoef[] = {
         0.00833333333333333333l,
        -0.00396825396825396825l,
         0.00416666666666666666l,
        -0.00757575757575757576l,
         0.02115090296908478727l,
        -0.08333333333333333333l
    };
    static constexpr int bounds[] = {
        6000, 320, 75, 30, 20, 10, 7
    };

    /// choose the amount of terms in series
    int degree = 0;
    while (x < bounds[degree])
        ++degree;

    /// apply
    long double firstTerm = taylorCoef[5] / x - 0.5;
    long double y = firstTerm / x;
    for (int i = 0; i < degree; ++i)
        y += taylorCoef[i] / std::pow(x, 2 * i + 4);
    return y;
}

double digamma(double x)
{
    /// Negative argument
    if (x < 0.0) {
        double y = 1.0 - x;
        double z = M_PI / std::tan(M_PI * y);
        return digamma(y) + z;
    }
    /// shift to minimum value,
    /// for which series expansion is applicable
    long double sum = 0;
    while (x < 7.0) {
        sum -= 1.0 / x;
        ++x;
    }
    return sum + digammamLogForLargeX(x) + std::log(x);
}

double digammamLog(double x)
{
    if (x < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    return (x < 7.0) ? digamma(x) - std::log(x) : digammamLogForLargeX(x);
}

double trigamma(double x)
{
    /// Negative argument
    if (x < 0.0)
    {
        double z = M_PI / std::sin(M_PI * x);
        return z * z - trigamma(1.0 - x);
    }
    /// set up tables
    static constexpr long double taylorCoef[] = {
         0.16666666666666666666l,
        -0.03333333333333333333l,
         0.02380952380952380952l,
        -0.03333333333333333333l,
         0.07575757575757575757l,
        -0.25311355311355311355l
    };
    static constexpr int bounds[] = {
        1000, 140, 47, 25, 15, 10
    };

    /// shift to minimum value,
    /// for which series expansion is applicable
    long double y = 0;
    while (x < 10.0) {
        y += 1.0 / (x * x);
        ++x;
    }

    /// choose the amount of terms in series
    int degree = 1;
    while (x < bounds[degree - 1])
        ++degree;

    /// apply
    long double firstTerm = 1.0 / x;
    firstTerm += 0.5 / (x * x);
    firstTerm += taylorCoef[0] / std::pow(x, 3);
    y += firstTerm;
    for (int i = 1; i < degree; ++i)
        y += taylorCoef[i] / std::pow(x, 2 * i + 3);
    return y;
}

enum REGULARISED_GAMMA_METHOD_ID {
    PT,
    QT,
    PUA,
    QUA,
    CF,
    UNDEFINED
};

double pgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId);
double lpgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId);
double qgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId);
double lqgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId);

REGULARISED_GAMMA_METHOD_ID getRegularizedGammaMethodId(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return UNDEFINED;
    double alpha = (x < 0.5) ? M_LN2 / (M_LN2 - logX) : x;
    if ((a <= 12 && a >= alpha) || 0.3 * a >= x)
        return PT;
    if (x <= 1.5 && a <= alpha)
        return QT;
    if (a <= 12 || 2.35 * a <= x)
        return CF;
    return (a > alpha) ? PUA : QUA;
}

double incompleteGammaUniformExpansion(double a, double x, double logX, double logA, bool isP)
{
    /// Uniform asymptotic expansion for P(a, x) if isP == true, or Q(a, x) otherwise
    static constexpr long double d[] = {-0.33333333333333333333l, 0.08333333333333333333l, -0.01481481481481481481l, 0.00115740740740740741l,
                                         0.00035273368606701940l, -0.000178755144032922l, 0.0000391926317852244l, -0.00000218544851067999l,
                                        -0.00000185406221071516l, 0.829671134095309e-6l, -0.176659527368261e-6l, 0.670785354340150e-8l,
                                         0.102618097842403e-7l, -0.438203601845335e-8l, 0.914769958223678e-9l, -0.255141939949460e-10l,
                                        -0.583077213255043e-10l, 0.243619480206674e-10l, -0.502766928011417e-11l, 0.110043920319559e-12l,
                                         0.337176326240099e-12l, -0.139238872241816e-12l, 0.285348938070474e-13l, -0.513911183424242e-15l,
                                        -0.197522882943494428e-16l, 0.809952115670456133e-17};
    static constexpr int N = 25;
    double lambda = x / a;
    double logLambda = logX - logA;
    double aux = x - a - a * logLambda;
    double eta = 0.0, base = 0.5;
    if (aux > 0.0) { /// otherwise, x ~ a and aux ~ 0.0
        eta = std::sqrt(2 * (lambda - 1.0 - logLambda));
        base = 0.5 * std::erfc(std::sqrt(aux));
    }
    long double sum = 0.0l;
    double betanp2 = d[N], betanp1 = d[N - 1];
    for (int n = N - 2; n >= 0; --n) {
        double beta = (n + 2) * betanp2 / a;
        beta += (x < a && (n & 1)) ? -d[n] : d[n];
        sum += beta * std::pow(eta, n);
        betanp2 = betanp1;
        betanp1 = beta;
    }
    sum *= a / (a + betanp2);
    double z = aux + 0.5 * (M_LN2 + M_LNPI + logA);
    double y = std::exp(-z) * sum;
    return base + (isP ? -y : y);
}

double lpgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == PT)
    {
        /// Taylor expansion of P(a,x)
        int n0 = 70.0 * x / a + 7; /// ~ from 7 to 28
        long double sum = 0.0;
        double lgammaAp1 = logA + lgammaA;
        for (int n = n0; n > 0; --n) {
            double addon = n * logX - std::lgamma(a + n + 1) + lgammaAp1;
            addon = std::exp(addon);
            sum += addon;
        }
        return a * logX - x + std::log1p(sum) - lgammaAp1;
    }
    return (mId == PUA) ? std::log(pgammaRaw(a, x, logX, logA, lgammaA, mId)) :
                          std::log1p(-qgammaRaw(a, x, logX, logA, lgammaA, mId));
}

double lpgamma(double a, double x, double logA, double lgammaA)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return RandMath::log1mexp(-x);
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lpgammaRaw(a, x, logX, logA, lgammaA, mId);
}

double lpgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return RandMath::log1mexp(-x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lpgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double lpgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return RandMath::log1mexp(-x);
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lpgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double pgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == PT)
        return std::exp(lpgammaRaw(a, x, logX, logA, lgammaA, mId));
    if (mId == PUA)
        return incompleteGammaUniformExpansion(a, x, logX, logA, true);
    return (mId == QUA) ? 1.0 - qgammaRaw(a, x, logX, logA, lgammaA, mId) :
                          -std::expm1(lqgammaRaw(a, x, logX, logA, lgammaA, mId));
}

double pgamma(double a, double x, double logA, double lgammaA)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 1.0;
    if (a == 1.0)
        return -std::expm1(-x);
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return pgammaRaw(a, x, logX, logA, lgammaA, mId);
}

double pgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 1.0;
    if (a == 1.0)
        return -std::expm1(-x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return pgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double pgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 1.0;
    if (a == 1.0)
        return -std::expm1(-x);
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return pgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double qtGammaExpansionAux(double a, double logX, double logA, double lgammaA)
{
    /// auxilary function for Taylor expansion of Q(a, x)
    long double sum = 0.0;
    for (int n = 1; n != 20; ++n) {
        double addon = n * logX - RandMath::lfact(n);
        addon = std::exp(addon);
        addon /= (a + n);
        sum += (n & 1) ? -addon : addon;
    }
    sum *= a;
    double y = std::log1p(sum);
    y += a * logX;
    y -= logA + lgammaA;
    return y;
}

double lqgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == QT)
    {
        double y = qtGammaExpansionAux(a, logX, logA, lgammaA);
        return RandMath::log1mexp(y);
    }
    if (mId == CF)
    {
        /// Continued fraction
        int k0 = std::min(40.0 / (x - 1) + 5.0, 60.0);
        long double sum = 0.0;
        double rhok = 0.0, tk = 1.0;
        for (int k = 1; k <= k0; ++k) {
            /// Calculate a(k)
            double ak = k * (a - k);
            double temp = x + 2 * k - a;
            ak /= temp * temp - 1;
            /// Calculate rho(k)
            ++rhok;
            rhok *= ak;
            rhok /= -(1 + rhok);
            /// Calculate t(k) and add it to the sum
            tk *= rhok;
            sum += tk;
        }
        double y = std::log1p(sum);
        y += a * logX;
        y -= x + lgammaA;
        y -= std::log1p(x - a);
        return y;
    }
    return (mId == QUA) ? std::log(qgammaRaw(a, x, logX, logA, lgammaA, mId)) :
                          std::log1p(-pgammaRaw(a, x, logX, logA, lgammaA, mId));
}

double lqgamma(double a, double x, double logA, double lgammaA)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    if (a == 1.0)
        return -x;
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lqgammaRaw(a, x, logX, logA, lgammaA, mId);
}

double lqgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    if (a == 1.0)
        return -x;
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lqgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double lqgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return -INFINITY;
    if (a == 1.0)
        return -x;
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return lqgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double qgammaRaw(double a, double x, double logX, double logA, double lgammaA, REGULARISED_GAMMA_METHOD_ID mId)
{
    if (mId == CF)
        return std::exp(lqgammaRaw(a, x, logX, logA, lgammaA, mId));
    if (mId == QT)
        return -std::expm1(qtGammaExpansionAux(a, logX, logA, lgammaA));
    if (mId == QUA)
        return incompleteGammaUniformExpansion(a, x, logX, logA, false);
    return (mId == PUA) ? 1.0 - pgammaRaw(a, x, logX, logA, lgammaA, mId) :
                          -std::expm1(lpgammaRaw(a, x, logX, logA, lgammaA, mId));
}

double qgamma(double a, double x, double logA, double lgammaA)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::exp(-x);
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return qgammaRaw(a, x, logX, logA, lgammaA, mId);
}

double qgamma(double a, double x, double logX)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::exp(-x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return qgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}

double qgamma(double a, double x)
{
    if (x < 0.0 || a < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;
    if (a == 1.0)
        return std::exp(-x);
    double logX = std::log(x);
    REGULARISED_GAMMA_METHOD_ID mId = getRegularizedGammaMethodId(a, x, logX);
    return qgammaRaw(a, x, logX, std::log(a), std::lgamma(a), mId);
}


}
