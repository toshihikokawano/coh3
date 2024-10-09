// coefficients of Gauss_Legendre integral

#ifdef GAUSS_LEGENDRE16
static const int MAX_GAUSSLEG = 8;
static double
gaussleg_x[] = {
  0.98940093499164993 , 0.94457502307323258 , 0.86563120238783174 ,
  0.75540440835500303 , 0.61787624440264375 , 0.45801677765722739 ,
  0.28160355077925891 , 0.095012509837637440},
gaussleg_a[] = {
  0.027152459411754095, 0.062253523938647893, 0.095158511682492785,
  0.12462897125553387 , 0.14959598881657673 , 0.16915651939500254 ,
  0.18260341504492359 , 0.18945061045506850 };
#endif


#ifdef GAUSS_LEGENDRE20
static const int MAX_GAUSSLEG = 10;
static double
gaussleg_x[] = {
  0.9931285991850949,  0.9639719272779138,  0.9122344282513259,
  0.8391169718222188,  0.7463319064601508,  0.6360536807265150,
  0.5108670019508271,  0.3737060887154196,  0.2277858511416451,
  0.0765265211334973},
gaussleg_a[] = {
  0.0176140071391521,  0.0406014298003869,  0.0626720483341091,
  0.0832767415767047,  0.1019301198172404,  0.1181945319615184,
  0.1316886384491766,  0.1420961093183821,  0.1491729864726037,
  0.1527533871307258};
#endif


#ifdef GAUSS_LEGENDRE30
static const int MAX_GAUSSLEG = 15;
static double
gaussleg_x[] = {
  0.996893484074649540272 ,  0.98366812327974720997  ,  0.960021864968307512217 ,
  0.926200047429274325879 ,  0.882560535792052681543 ,  0.829565762382768397443 ,
  0.767777432104826194918 ,  0.697850494793315796932 ,  0.6205261829892428611405,
  0.5366241481420198992642,  0.4470337695380891767806,  0.352704725530878113471 ,
  0.2546369261678898464398,  0.1538699136085835469638,  0.051471842555317695833},
gaussleg_a[] = {
  0.007968192496166605615 ,  0.018466468311090959142 ,  0.0287847078833233693497,
  0.038799192569627049597 ,  0.048402672830594052903 ,  0.057493156217619066482 ,
  0.065974229882180495128 ,  0.073755974737705206268 ,  0.0807558952294202153547,
  0.0868997872010829798024,  0.0921225222377861287176,  0.0963687371746442596395,
  0.099593420586795267063 ,  0.101762389748405504596 ,  0.102852652893558840341};
#endif


#ifdef GAUSS_LEGENDRE32
static const int MAX_GAUSSLEG = 16;
static double
gaussleg_x[] = {
  0.997263861849481563545 ,  0.9856115115452683354002,  0.9647622555875064307738,
  0.9349060759377396891709,  0.8963211557660521239653,  0.8493676137325699701337,
  0.7944837959679424069631,  0.7321821187402896803874,  0.6630442669302152009751,
  0.5877157572407623290407,  0.5068999089322293900237,  0.4213512761306353453641,
  0.3318686022821276497799,  0.2392873622521370745446,  0.1444719615827964934852,
  0.0483076656877383162348},
gaussleg_a[] = {
 0.0070186100094700966004 , 0.0162743947309056706052 , 0.0253920653092620594558 ,
 0.0342738629130214331027 , 0.04283589802222668065688, 0.05099805926237617619616,
 0.0586840934785355471453 , 0.06582222277636184683765, 0.0723457941088485062254 ,
 0.0781938957870703064717 , 0.0833119242269467552222 , 0.087652093004403811143  ,
 0.0911738786957638847129 , 0.093844399080804565639  , 0.0956387200792748594191 ,
 0.0965400885147278005668};
#endif


#ifdef GAUSS_LEGENDRE40
static const int MAX_GAUSSLEG = 20;
static double
gaussleg_x[] = {
  0.9982377097105592003496,  0.9907262386994570064531,  0.9772599499837742626634,
  0.9579168192137916558045,  0.9328128082786765333609,  0.9020988069688742967283,
  0.8659595032122595038208,  0.8246122308333116631963,  0.778305651426519387695 ,
  0.727318255189927103281 ,  0.6719566846141795483794,  0.6125538896679802379526,
  0.5494671250951282020759,  0.4830758016861787129086,  0.4137792043716050015249,
  0.3419940908257584730075,  0.2681521850072536811412,  0.1926975807013710997155,
  0.1160840706752552084835,  0.0387724175060508219332},
gaussleg_a[] = {
  0.0045212770985331912585,  0.01049828453115281361474, 0.0164210583819078887129,
  0.0222458491941669572615,  0.0279370069800234010985,  0.0334601952825478473927,
  0.03878216797447201764  ,  0.0438709081856732719917,  0.0486958076350722320614,
  0.053227846983936824355 ,  0.057439769099391551367 ,  0.0613062424929289391665,
  0.0648040134566010380746,  0.0679120458152339038257,  0.0706116473912867796955,
  0.0728865823958040590605,  0.0747231690579682642002,  0.0761103619006262423716,
  0.0770398181642479655883,  0.077505947978424811264};
#endif


#ifdef GAUSS_LEGENDRE50
static const int MAX_GAUSSLEG = 25;
static double
gaussleg_x[]={
  0.9988664044200710501855 ,  0.99403196943209071258518,  0.985354084048005882309  ,
  0.9728643851066920737133 ,  0.9566109552428079429977 ,  0.9366566189448779337809 ,
  0.91307855665579189308979,  0.8859679795236130486375 ,  0.8554297694299460846114 ,
  0.8215820708593359483563 ,  0.7845558329003992639053 ,  0.74449430222606853826051,
  0.7015524687068222510895 ,  0.65589646568543936078163,  0.6077029271849502391804 ,
  0.5571583045146500543155 ,  0.5044581449074642016515 ,  0.4498063349740387891471 ,
  0.3934143118975651273942 ,  0.335500245419437356837  ,  0.2762881937795319903276 ,
  0.21600723687604175684735,  0.1548905899981459020716 ,  0.09317470156008614085445,
  0.03109833832718887611233},
gaussleg_a[]={
  0.0029086225531551409584 ,  0.00675979919574540150278,  0.01059054838365096926357,
  0.0143808227614855744194 ,  0.0181155607134893903513 ,  0.0217802431701247929816 ,
  0.02536067357001239044019,  0.0288429935805351980299 ,  0.0322137282235780166482 ,
  0.035459835615146154161  ,  0.0385687566125876752448 ,   0.04152846309014769742241,
  0.044327504338803275492  ,  0.04695505130394843296563,  0.0494009384494663149212 ,
  0.0516557030695811384899 ,  0.0537106218889962465235 ,  0.0555577448062125176236 ,
  0.057189925647728383723  ,  0.0586008498132224458351 ,  0.0597850587042654575096 ,
  0.06073797084177021603175,  0.0614558995903166637564 ,  0.0619360674206832433841 ,
  0.06217661665534726232103};
#endif


#ifdef GAUSS_LEGENDRE60
static const int MAX_GAUSSLEG = 30;
static double
gaussleg_x[]={
  0.9992101232274360220342,  0.9958405251188381738767,  0.9897878952222217173673,
  0.9810672017525981856186,  0.9697017887650527337215,  0.9557222558399961073972,
  0.9391662761164232494954,  0.9200784761776275528567,  0.8985103108100459419378,
  0.8745199226468983151293,  0.8481719847859296324905,  0.8195375261621457593685,
  0.7886937399322640545699,  0.7557237753065856868688,  0.720716513355730399436 ,
  0.6837663273813554372229,  0.6449728284894770678134,  0.6044405970485103634442,
  0.5622789007539445391783,  0.5186014000585697474179,  0.4735258417617071111082,
  0.4271737415830783893075,  0.379670056576797977155 ,  0.3311428482684481942524,
  0.2817229374232616916907,  0.2315435513760293380103,  0.1807399648734254172409,
  0.1294491353969450031464,  0.07780933394953656941929,  0.02595977230124779858917},
gaussleg_a[]={
  0.0020268119688737584964 ,  0.00471272992695356864089,  0.0073899311633454555315 ,
  0.01004755718228798435789,  0.0126781664768159601315 ,  0.0152746185967847993067 ,
  0.0178299010142077202604 ,  0.020337120729457286775  ,  0.0227895169439978198638 ,
  0.0251804776215212483796 ,  0.0275035567499247916352 ,  0.0297524915007889452408 ,
  0.0319212190192963289495 ,  0.0340038927249464228349 ,  0.0359948980510845030666 ,
  0.03788886756924344403094,  0.0396806954523807994701 ,  0.0413655512355847556132 ,
  0.0429388928359356419542 ,  0.0443964787957871133278 ,  0.0457343797161144866472 ,
  0.046948988848912204847  ,  0.04803703181997118096367,  0.04899557545575683538948,
  0.0498220356905501810112 ,  0.0505141845325093745982 ,  0.0510701560698556274045 ,
  0.051488451500980933995  ,  0.0517679431749101875438 ,  0.0519078776312206397329};
#endif





#ifdef GAUSS_LAGUERRE32
static const int MAX_GAUSSLAG = 32;
static double
gausslag_x[] = {
  0.0444893658332670184189,  0.234526109519618537453 ,  0.576884629301886426492 ,
  1.072448753817817633041 ,  1.722408776444645441131 ,  2.528336706425794881124 ,
  3.492213273021994489609 ,  4.616456769749767387762 ,  5.903958504174243946562 ,
  7.358126733186241113222 ,  8.982940924212596103378 ,  10.7830186325399720675  ,
  12.76369798674272511497 ,  14.9311397555225573198  ,  17.29245433671531478924 ,
  19.85586094033605473979 ,  22.63088901319677448868 ,  25.62863602245924776748 ,
  28.86210181632347474434 ,  32.34662915396473700323 ,  36.10049480575197380402 ,
  40.14571977153944153621 ,  44.50920799575493797591 ,  49.22439498730863917672 ,
  54.33372133339690733287 ,  59.89250916213401819613 ,  65.97537728793505279656 ,
  72.68762809066270863868 ,  80.18744697791352306749 ,  88.73534041789239868936 ,
  98.82954286828397255918 ,  111.7513980979376952137 },

gausslag_a[] = {
 0.1092183419523849711      , 0.2104431079388132329      , 0.2352132296698480054      ,
 0.19590333597288104341     , 0.12998378628607176061     , 0.07057862386571744156     ,
 0.031760912509175070306    , 0.011918214834838557057    , 0.0037388162946115247897   ,
 9.808033066149551322e-04   , 2.14864918801364188e-04    , 3.9203419679879472043e-05  ,
 5.93454161286863287836e-06 , 7.4164045786675522191e-07  , 7.6045678791207814811e-08  ,
 6.35060222662580674243e-09 , 4.2813829710409288788e-10  , 2.30589949189133607927e-11 ,
 9.7993792887270940633e-13  , 3.23780165772926646231e-14 , 8.1718234434207194332e-16  ,
 1.54213383339382337218e-17 , 2.11979229016361861204e-19 , 2.054429673788045426656e-21,
 1.34698258663739515581e-23 , 5.66129413039735937113e-26 , 1.41856054546303690595e-28 ,
 1.913375494454224309371e-31, 1.192248760098222356542e-34, 2.671511219240136986e-38   ,
 1.33861694210625628272e-42 , 4.51053619389897423222e-48 };
#endif


#ifdef GAUSS_LAGUERRE40
static const int MAX_GAUSSLAG = 40;
static double
gausslag_x[] = {
 0.035700394308888385122  , 0.188162283158698516004  , 0.4626942813145764535649 ,
 0.8597729639729349222573 , 1.380010820527337186498  , 2.024209135922826733442  ,
 2.793369353506816457654  , 3.688702677908270209592  , 4.711641146554972693619  ,
 5.863850878343718114273  , 7.147247908102288250686  , 8.564017017586163762719  ,
 10.11663404845193940685  , 11.80789229400458484284  , 13.64093371253708722837  ,
 15.6192858933390738372   , 17.74690595009566304257  , 20.02823283457489052961  ,
 22.46824998349841835137  , 25.0725607724262037944   , 27.84748000916886272075  ,
 30.80014573944546270075  , 33.9386570849137196091   , 37.27224588047600432832  ,
 40.81149282388692046616  , 44.56860317533446270712  , 48.55776353305999228096  ,
 52.79561118721693296935  , 57.30186332339362749503  , 62.10017907277511161217  ,
 67.21937092712699879908  , 72.69515884761246211752  , 78.57280291157130928054  ,
 84.9112311357049845427   , 91.78987467123637699234  , 99.32080871744680825011  ,
 107.6724406393882725208  , 117.1223095126906888077  , 128.2018419882556511925  ,
 142.2800444691599978883  },

gausslag_a[] = {
 0.0884121061903424409       , 0.1768147390957222956       , 0.211363117015962431        ,
 0.19408119531860179966      , 0.146434282424125114407     , 0.093326798435770880507     ,
 0.05093220436104423703      , 0.02397619301568484184      , 0.009774625246714459619     ,
 0.00345793999301848686132   , 0.001062246893896871935     , 2.8327168532432471583E-4    ,
 6.55094050032462928E-5      , 1.3116069073267784125E-5    , 2.26845287877936505455E-6   ,
 3.37962648220067921081E-7   , 4.3228213222820885689E-8    , 4.7284937709907793279E-9    ,
 4.4031741042328488129E-10   , 3.47244148480382248556E-11  , 2.30538154491682216159E-12  ,
 1.279772597676635607209E-13 , 5.8941771723511529447E-15   , 2.23221757990457741835E-16  ,
 6.88033648428430234095E-18  , 1.70560373681808674846E-19  , 3.353711940666182935466E-21 ,
 5.14619956013667914084E-23  , 6.04476251158766328897E-25  , 5.31058477732133075276E-27  ,
 3.39252805328052189606E-29  , 1.521735493181456997543E-31 , 4.58529161450268691756E-34  ,
 8.762158657486248561E-37    , 9.8274157251479333061E-40   , 5.80115201916977910848E-43  ,
 1.53090868460668685357E-46  , 1.3819863056493280997E-50   , 2.566633605012372183814E-55 ,
 2.7003609402170336406E-61   };
#endif



