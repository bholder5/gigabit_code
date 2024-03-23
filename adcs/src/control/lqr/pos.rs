extern crate nalgebra as na;

pub type AMat1 = na::SMatrix<f64, 18, 18>;

pub fn init_AMat1_pos() -> AMat1 {
    let a = AMat1::from_row_slice(&[
0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0689462086337438,-0.0157519186630494,0.0116656804601277,-0.3714767823551953,-0.0234946288222179,0.0254953098606162,-0.0004906064405125,-0.0029675000579011,-0.0002122127468680,-0.0020861960680826,0.0007938485810442,0.0024033891416388,-0.0001444222476027,-0.0003196796102788,-0.0002172718102152,-0.0004798412018590,0.0002762971705204,0.0010847231401333,
-0.0027754168497183,-1.6143807506811425,0.0001452228171314,-0.0149537290223547,-2.4079147008269470,0.0003173840338120,0.0085573602661202,0.0126736647787327,0.0248494430201144,0.0369927351720856,0.0496872316103949,0.0741597041682725,-0.0001737565710798,-0.0002606392541449,-0.0001601025919078,-0.0002409437119157,0.0104644549362131,0.0156351987124577,
0.0077792750246400,0.0005496283764197,-0.8339603178886783,0.0419141256999456,0.0008197931293557,-1.8226177880237457,0.0001078123678953,0.0004449420695222,0.0001010318877689,0.0003907392708709,0.0002763881115069,0.0005024103760820,0.0100885312354158,0.0220490020938188,0.0152876303098661,0.0334117170160411,-0.0001652340684415,-0.0004208950944746,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.5597552222732161,2.0292953163990823,0.0055649123269603,-3.0159173796077745,3.0267766898329000,0.0121620992969889,-6850.1034513899094236,-0.2067132741338955,-0.0349456262187941,-0.0664859826903999,-0.0598964969454370,-0.0794142672265354,0.0001355988078258,0.0001040768405018,0.0000804353551190,-0.0000106775774847,-0.0117456222204423,-0.0120688529104336,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.4430140280308673,5.8812648700206536,0.0070484872560834,-2.3869249510926265,8.7721463068763263,0.0154044478808297,-0.0349456262187941,-0.0664859826903999,-7962.7641700026651961,-0.3293126329115477,-0.1789538359389538,-0.2590635486549837,0.0005306480441107,0.0006915868314917,0.0004319207662061,0.0005071888722943,-0.0369890161790715,-0.0508567760951567,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.3061401830938242,11.7481498981390935,0.0397615351110378,1.6494593744733512,17.5228444933527427,0.0868987164233178,-0.0598964969454370,-0.0794142672265354,-0.1789538359389538,-0.2590635486549837,-12648.2293973437281238,-0.7716457237244130,0.0007774322185984,0.0008539300746463,0.0004247013457750,0.0001603932698407,-0.0768592948894681,-0.1176145292118073,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0015872087984565,-0.0406627962594211,1.2203331397338266,-0.0085517569284863,-0.0606502182638635,2.6670344380709956,0.0001355988078258,0.0001040768405018,0.0005306480441107,0.0006915868314917,0.0007774322185984,0.0008539300746463,-13102.1494961056614557,-0.2611984780865440,-0.0223739876328379,-0.0488956721455619,0.0004755018374104,0.0008691187637550,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0018463204278117,-0.0372169598996067,1.8492408055749974,-0.0099478301318004,-0.0555106128616433,4.0415102664756040,0.0000804353551190,-0.0000106775774847,0.0004319207662061,0.0005071888722943,0.0004247013457750,0.0001603932698407,-0.0223739876328379,-0.0488956721455619,-17096.5677736761390406,-0.3355982153907916,0.0005610915411643,0.0010737380280361,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,
0.1680882682356104,2.4732172777472581,-0.0172056473442929,0.9056464492126054,3.6889043919249880,-0.0376028909667592,-0.0117456222204423,-0.0120688529104336,-0.0369890161790715,-0.0508567760951567,-0.0768592948894681,-0.1176145292118073,0.0004755018374104,0.0008691187637550,0.0005610915411643,0.0010737380280361,-19067.2922443990864849,-0.3023335676739867,
    ]);
    return a;
}

pub type AMat2 = na::SMatrix<f64, 18, 5>;

pub fn init_AMat2_pos() -> AMat2 {
    let a = AMat2::from_row_slice(&[
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0010248848411720,0.0000412565493769,0.0000412565493769,-0.0001156388613135,-0.0001156388613135,
0.0000412565493769,0.0042282962843103,0.0042282962843103,-0.0000014395560779,-0.0000014395560779,
-0.0001156388613135,-0.0000014395560779,-0.0000014395560779,0.0082668320866901,0.0082668320866901,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0083207569230769,-0.0202840893682008,0.0096540544864420,-0.0040428003580101,0.0039324732961646,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0065854000000000,-0.0191757382236295,-0.0116320255414389,0.0074533822429512,-0.0075931218908971,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0045507713846154,-0.0503321172932141,-0.0112080844724667,0.0067533698625509,-0.0075416614088556,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000235938461538,0.0001740123658819,0.0000389911034182,-0.0121435815217691,-0.0120501095679042,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000274455384615,0.0000075046009171,0.0001874485824510,-0.0184522418228579,-0.0182098460478841,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0024986307692308,-0.0331079744840978,0.0201525481742034,-0.0053156546824437,0.0056567649115214,
    ]);
    return a;
}

pub type AMat3 = na::SMatrix<f64, 18, 18>;

pub fn init_AMat3_pos() -> AMat3 {
    let a = AMat3::from_row_slice(&[
0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0689462160491024,-0.0157519348100180,0.0116656814374637,-0.3714768066196227,-0.0234946224459146,0.0254953103936908,-0.0004782276041792,-0.0029870234384987,-0.0002596187110628,-0.0021193550988824,0.0007672845054729,0.0024034387949698,-0.0001444441864434,-0.0003197337889087,-0.0002173826961869,-0.0004799048987665,0.0002731232407649,0.0010847751288962,
-0.0027753935232604,-1.6143807624071471,0.0001452249206424,-0.0149535499131652,-2.4079153267704951,0.0003173878220022,0.0083032685703718,0.0127810130344782,0.0262124614852005,0.0377772379409067,0.0502588432752922,0.0741613306143083,-0.0001753526839825,-0.0002606564981992,-0.0001608732675516,-0.0002409566455387,0.0105394637144998,0.0156361002687906,
0.0077792779138375,0.0005496146984148,-0.8339603998149754,0.0419141367855421,0.0008197630496848,-1.8226178713127235,0.0001185296148191,0.0004494182464672,0.0001354290480931,0.0004020302699284,0.0002852534697821,0.0005024472640849,0.0100938323702569,0.0220528300764659,0.0152972757475928,0.0334162317174673,-0.0001643490871968,-0.0004209469429236,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.5597552573653369,2.0292943635816232,0.0055649007700793,-3.0159176242302252,3.0267769132189541,0.0121620735250808,-6850.1026918851457594,-0.2081729711472061,-0.0374994286523080,-0.0681234843224870,-0.0612967251208331,-0.0794195068059423,0.0001392113068345,0.0001040989188547,0.0000816513544473,-0.0000107507117947,-0.0119258830406481,-0.0120835496372758,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.4430140562016650,5.8812632553827378,0.0070484747859780,-2.3869252559378560,8.7721470248413063,0.0154044218798121,-0.0333691653210492,-0.0671393900563515,-7962.7690413426726082,-0.3332935153861976,-0.1818671205366786,-0.2590709459292178,0.0005385348584590,0.0006916345875211,0.0004356374092948,0.0005072091655935,-0.0373672404527156,-0.0508617410191321,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.3061402529172211,11.7481470182966614,0.0397615108913762,1.6494591097877558,17.5228459964098953,0.0868986681489150,-0.0570114752233489,-0.0813841780350032,-0.1882394776948886,-0.2661431723667216,-12648.2347179342814343,-0.7716634373064212,0.0007919565211768,0.0008539320598088,0.0004309021951225,0.0001602412225369,-0.0775548946805199,-0.1176372054570564,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0015872123182935,-0.0406627602712297,1.2203332597484478,-0.0085517694515690,-0.0606501772493753,2.6670345602434735,0.0001064772745587,0.0001070883627715,0.0005197763311282,0.0007075790869740,0.0007889574726986,0.0008539332467949,-13102.1495039274377632,-0.2612040801635917,-0.0223881314302648,-0.0489022787283596,0.0004773594437892,0.0008692618211946,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0018463253225445,-0.0372169159218156,1.8492409873077893,-0.0099478491062832,-0.0555105499930591,4.0415104513511748,0.0000441568235249,-0.0000224262061491,0.0003908915896992,0.0005157016254484,0.0004283259235500,0.0001603285094311,-0.0223858027019262,-0.0489041606452061,-17096.5677950968674850,-0.3356082272386405,0.0005621004903945,0.0010737623284548,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,
0.1680883685790757,2.4732155960874445,-0.0172056663284828,0.9056466670450841,3.6889038651038901,-0.0376029263832161,-0.0107382402480112,-0.0142115804947686,-0.0394421026538164,-0.0526745286923377,-0.0785021220670676,-0.1176234399562198,0.0004800933196982,0.0008692815928321,0.0005626386556114,0.0010737591197905,-19067.2924586703156820,-0.3023595083741275,
    ]);
    return a;
}

pub type AMat4 = na::SMatrix<f64, 18, 5>;

pub fn init_AMat4_pos() -> AMat4 {
    let a = AMat4::from_row_slice(&[
0.0000001649104306,-0.0000004617748504,-0.0000002044685653,0.0000001385313385,-0.0000001410844145,
-0.0000002192953786,0.0000006418897641,0.0000003473375861,-0.0000002278642820,0.0000002323694394,
0.0000000027223436,-0.0000000075599340,-0.0000000032481649,0.0000000024899368,-0.0000000018860470,
0.0000006163922399,0.0000012074618795,0.0000010860103379,-0.0000002527069384,0.0000000418895085,
-0.0000045929775675,0.0000176039981087,0.0000051879789102,-0.0000098161058755,0.0000097565237676,
-0.0000001162095341,0.0000000048052479,-0.0000000806179511,0.0000081864303039,0.0000082845064026,
0.0000334430642489,-0.0000943756571878,-0.0000675179314510,0.0000306840438972,-0.0000324906965680,
0.0835422301871644,-0.2036709145759154,0.0968728467780380,-0.0404510848839790,0.0393466573600329,
0.0000770642713730,-0.0002235993234913,-0.0002907456051264,-0.0000491719813324,0.0000437713610965,
0.5297662519867699,-1.5422841082241658,-0.9358190437780041,0.5998440894072899,-0.6110725536515889,
-0.0000312437946938,0.0000877758555100,0.0000510706915543,-0.0000335804227324,0.0000342610525091,
-0.0007608564580149,-0.0106938447463579,-0.0026190155756841,0.0014874100748725,-0.0016519080121183,
0.0000001045800835,-0.0000001216007128,-0.0000000946111004,-0.0000016077695927,-0.0000017544216437,
0.0000231381805968,0.0001722374006107,0.0000389052504397,-0.0121403318508314,-0.0120473025847040,
0.0000000589110679,-0.0000000521995081,-0.0000000490415152,-0.0000019157017642,-0.0000019890623532,
0.0000273521428561,0.0000076546179239,0.0001875341355721,-0.0184487519577985,-0.0182063803428103,
-0.0000027995987010,0.0000040673790648,0.0000022418113127,-0.0000024675995620,0.0000026066087229,
-0.0049894522602494,-0.0662304712730999,0.0402977136551287,-0.0106305540607511,0.0113123507155264,
    ]);
    return a;
}

pub type AMat5 = na::SMatrix<f64, 5, 18>;

pub fn init_AMat5_pos() -> AMat5 {
    let a = AMat5::from_row_slice(&[
0.0000074324948759,0.0000158612506227,0.0000001848223955,0.0000252781508232,-0.0000117905079576,0.0000006538484642,-0.0146713000485754,0.0200492586891237,0.0588715840201857,0.0397432304824944,0.0312990602033500,-0.0000335404648951,-0.0000662570787314,0.0000004595124550,-0.0000308318796937,0.0000004159665646,0.0038060485458172,-0.0000415217621658,
0.0000005388655344,-0.0000475102322057,-0.0000005488917421,-0.0000075488100552,0.0000289127112278,-0.0000010577428097,0.0480920990602725,-0.0488690148868669,-0.1640331776549474,-0.1157054897974260,-0.0897946417930220,-0.0003210960784763,0.0002471618609678,0.0000032169344142,0.0001050468785976,0.0000000852708568,-0.0117216182116667,-0.0005469520779498,
-0.0000061282212170,0.0000501293310841,0.0000000529789487,-0.0000350578330853,0.0001192403460084,0.0000001588821422,0.0121437073355384,0.0232852350732476,-0.1588988254361222,-0.0702189575416086,-0.0456982056986110,-0.0000632358978177,0.0001307497076251,0.0000006991856070,0.0000771228210393,0.0000027835671627,-0.0060552544538010,0.0003341393778640,
-0.0000055559301201,0.0000795200386830,0.0000050917120281,-0.0000280279530989,0.0000882385323062,0.0000053055407368,-0.0284258456620430,-0.0097621146151052,-0.0356158009123262,0.0450094691756476,0.0298228027436082,0.0000456839153360,-0.0004031330780808,-0.0002324079445573,-0.0006290210591755,-0.0002748641276480,0.0038229650180594,-0.0000877412637009,
0.0000053094320224,-0.0000776431465419,0.0000048210278231,0.0000270331591372,-0.0000847390654827,0.0000047785269040,0.0269346932586991,0.0094966522514419,0.0322222160641056,-0.0458517247984512,-0.0304809774940013,-0.0000506821844505,-0.0002389814167630,-0.0002306380967702,-0.0005381419780660,-0.0002712518629813,-0.0038798725401831,0.0000933952510321,
    ]);
    return a;
}

