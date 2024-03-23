extern crate nalgebra as na;

pub type AMat1 = na::SMatrix<f64, 18, 18>;

pub fn init_AMat1_pos() -> AMat1 {
    let a = AMat1::from_row_slice(&[
0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.1537327263726129,-0.0087714588860896,-0.0245857090076900,-0.3074654527452259,-0.0123769648130835,-0.0346916583939606,-0.0012311959557878,-0.0024906477963196,-0.0008755591804748,-0.0018317057883981,0.0009791915201554,0.0017937190105524,0.0002928482025947,0.0004110871724740,0.0004456074226314,0.0006262892949563,0.0004367470965510,0.0008425002287218,
-0.0061884824065417,-0.8989682262841721,-0.0003060606653399,-0.0123769648130835,-1.2684888852918752,-0.0004318668233843,0.0047265558270562,0.0066390786565318,0.0138068681261974,0.0194581758687224,0.0276896260143995,0.0390880259249756,-0.0000921853170004,-0.0001301640606131,-0.0000821877934591,-0.0001160710958774,0.0058386688376245,0.0082477589471756,
-0.0173458291969803,-0.0003060606653399,-1.7575919184183779,-0.0346916583939606,-0.0004318668233843,-2.4800496260117830,-0.0000457487344617,-0.0001495571828216,0.0000132879367002,-0.0000485253908375,0.0007811018462722,0.0011486634693282,0.0212608761171430,0.0299999127651776,0.0322179887776074,0.0454609044281717,-0.0002544428952109,-0.0003335059935246,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-1.2481135384615389,1.1300134806614612,0.0117281898464117,-2.4962270769230779,1.5945052322638322,0.0165490592768331,-6850.1043996934813549,-0.1947767516731249,-0.0256267405245252,-0.0410013803311026,-0.0290953720397957,-0.0377098508455750,-0.0000509735040332,-0.0000892694384575,-0.0001390952488293,-0.0002164448934776,-0.0041993279954923,-0.0040887752105801,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.9878099999999993,3.2749834549743904,0.0148548605642659,-1.9756199999999986,4.6211645647602602,0.0209609471918925,-0.0256267405245252,-0.0410013803311026,-7962.7276113964189790,-0.2626631810661392,-0.0962821287982605,-0.1332113665597436,0.0001457882278998,0.0001919882272581,0.0000198181460074,0.0000119973727869,-0.0187436973759232,-0.0249946585752665,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.6826157076923076,6.5419594922664128,0.0837984149557547,1.3652314153846152,9.2310302648521230,0.1182437319457048,-0.0290953720397957,-0.0377098508455750,-0.0962821287982605,-0.1332113665597436,-12648.0709329507171788,-0.5152246614541263,-0.0003008609186956,-0.0004150443661993,-0.0008796897563675,-0.0012302521661094,-0.0440682494247250,-0.0631870083287032,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0035390769230768,-0.0226430857860906,2.5718821605377693,-0.0070781538461536,-0.0319505203950189,3.6290536634509918,-0.0000509735040332,-0.0000892694384575,0.0001457882278998,0.0001919882272581,-0.0003008609186956,-0.0004150443661993,-13102.1658437123332988,-0.2728327847315578,-0.0471475891841142,-0.0665276508401135,0.0005941659176307,0.0008436056624970,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0041168307692270,-0.0207242711575460,3.8973205623459442,-0.0082336615384540,-0.0292429775052049,5.4993131806113027,-0.0001390952488293,-0.0002164448934776,0.0000198181460074,0.0000119973727869,-0.0008796897563675,-0.0012302521661094,-0.0471475891841142,-0.0665276508401135,-17096.6053155708068516,-0.3623186571330788,0.0008092402788411,0.0011479364895428,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,
0.3747946153846155,1.3772115087805299,-0.0362613257187646,0.7495892307692310,1.9433139464841560,-0.0511665343616508,-0.0041993279954923,-0.0040887752105801,-0.0187436973759232,-0.0249946585752665,-0.0440682494247250,-0.0631870083287032,0.0005941659176307,0.0008436056624970,0.0008092402788411,0.0011479364895428,-19067.2856645209794806,-0.2906385174633941,
    ]);
    return a;
}

pub type AMat2 = na::SMatrix<f64, 18, 5>;

pub fn init_AMat2_pos() -> AMat2 {
    let a = AMat2::from_row_slice(&[
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.0010248848411720,0.0000412565493769,0.0000412565493769,0.0001156388613135,0.0001156388613135,
0.0000412565493769,0.0042282962843103,0.0042282962843103,0.0000014395560779,0.0000014395560779,
0.0001156388613135,0.0000014395560779,0.0000014395560779,0.0082668320866901,0.0082668320866901,
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
-0.1699715926668470,-0.0094384542864944,-0.0253119839717273,-0.3570057768862296,-0.0126853051095947,-0.0347047924830508,-0.0012559706125616,-0.0025113719572100,-0.0013022094048650,-0.0022241433026011,0.0009377558471519,0.0018005772709813,0.0003247555873054,0.0004161453467301,0.0004955505755747,0.0006320783479201,0.0004693154931691,0.0008476516171639,
-0.0069088074989688,-0.9685908748198502,-0.0003127177175080,-0.0121389603688342,-1.3466062158769929,-0.0004337756406347,0.0031408313411340,0.0067419204153213,0.0234585674454631,0.0277711180825008,0.0340997781517743,0.0393165948383676,-0.0000999841831751,-0.0001316157258690,-0.0000891215134151,-0.0001170336772962,0.0070918128590857,0.0083331512230799,
-0.0201084091001131,-0.0003288824078760,-1.8149707453940591,-0.0351378049653520,-0.0004537863379289,-2.5302390628304541,-0.0000707055251069,-0.0001496551475776,-0.0001643933971426,-0.0000179302692625,0.0008628085854444,0.0011555449356456,0.0234798166490964,0.0303657606972816,0.0357399373526219,0.0458782571091693,-0.0002873922773097,-0.0003378055761399,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-1.3797490575779698,1.2176200395405734,0.0127316665353345,-2.9020350564195647,1.6964073952209149,0.0225145274777463,-6850.1015999884293706,-0.1962074576336118,-0.0427791094541459,-0.0583275170307463,-0.0389843804085565,-0.0384205088505284,-0.0000644982701361,-0.0000874614363618,-0.0001750962127273,-0.0002229616990267,-0.0059943568323417,-0.0054768339348828,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-1.0918109374411449,3.5286927841854645,0.0158303505769511,-2.3027283610693141,4.9087314590563462,0.0258971511074013,-0.0172749770799688,-0.0416269069382847,-7962.7630457404729896,-0.3058545054893864,-0.1239002144687668,-0.1342481992492643,0.0001588123335565,0.0001962733680553,0.0000179850707664,0.0000136491942308,-0.0237338720522291,-0.0254483052731338,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
0.7551965850538460,7.0485401845705544,0.0861930838140948,1.5691161291022278,9.7976042580312850,0.1177280027995954,-0.0128723809809521,-0.0395823173156658,-0.1605434342649500,-0.2087431492919232,-12648.1232399594991875,-0.5176797100825455,-0.0003730235294834,-0.0004156059994184,-0.0010288222788102,-0.0012444838764778,-0.0544547679888435,-0.0653146252093105,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0025507245347212,-0.0243956059303784,2.6558589232739385,-0.0157054379830864,-0.0338398288557200,3.7026251179713685,-0.0000791375988251,-0.0000868205022398,0.0005347408425745,0.0003539460294783,-0.0002458363991983,-0.0004154218018765,-13102.1690906418225495,-0.2733682002267331,-0.0523011981983806,-0.0671383967855425,0.0006837506675129,0.0008576326552626,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,0.0000000000000000,0.0000000000000000,
-0.0024865690802065,-0.0223278357053665,4.0245747438479089,-0.0209385041723798,-0.0309286961962064,5.6107949471543703,-0.0001462006725143,-0.0002285399754514,0.0004575579268161,0.0000959433331426,-0.0009052049685482,-0.0012402385349705,-0.0520679464694399,-0.0673389384588885,-17096.6131253828134504,-0.3632441998339925,0.0009208351577319,0.0011509732683943,
0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000000,
0.4144378336611866,1.4838228956392041,-0.0376321685126091,0.8672825614396455,2.0618938943053879,-0.0538542999627239,0.0002704810530770,-0.0061603680662530,-0.0330544409071906,-0.0436337059976795,-0.0559905405405348,-0.0643891439329798,0.0006203844709748,0.0008586790366671,0.0008421159087004,0.0011503758292256,-19067.2884165388750262,-0.2930918809503894,
    ]);
    return a;
}

pub type AMat4 = na::SMatrix<f64, 18, 5>;

pub fn init_AMat4_pos() -> AMat4 {
    let a = AMat4::from_row_slice(&[
-1.1640484957596984,-0.0583138120921469,-0.0583683206454149,-0.1086232795629266,-0.1086305200248291,
-0.0427833491364733,-7.8636227333646769,-7.8632073730774916,-0.0006648767240653,-0.0008080535115366,
-0.1976748892508887,-0.0020325957641834,-0.0020412282926233,-17.7892763433874954,-17.7892731430464686,
11.7521883058882306,0.4943276928229681,0.4943295197977279,1.1925001030888924,1.1924928301792046,
0.5109808396937348,49.2380501552573833,49.2396747582547434,0.0178328981001574,0.0172399829846069,
1.4209461022038472,0.0251900375926496,0.0252329170575065,85.2975493661585062,85.2975338176630373,
-0.0043237484805562,0.0098047220339766,0.0101180747675994,0.0001016236492903,-0.0005709111047713,
0.0835454146955288,-0.1928061384087535,0.0815334905266923,-0.0388061664610609,0.0387390150200309,
-0.0023679788785545,0.0253823109948996,0.0253140123302443,-0.0003523346884920,0.0002222218542690,
5.6670513900185382,-15.0161365410533012,-9.1167476743665290,6.4161332310229184,-6.5360472431899526,
0.0015168375406619,0.0343619238601313,0.0341956686498352,0.0005327155536607,0.0011212250352854,
-0.0906307475155817,-0.9378223910309725,-0.2236624280083393,0.1345708722318865,-0.1533105873862109,
0.0003432454676637,-0.0001092189638704,-0.0001070404672083,0.0202645765713630,0.0202626031963232,
0.0016474363533522,0.0159676105260624,0.0036101372775846,-1.2492299496006352,-1.2398427944411301,
0.0004010669336577,-0.0000704272033432,-0.0000723714056802,0.0235345007919426,0.0235334160734875,
0.0019120459447950,0.0007146642533024,0.0170762071721657,-1.8866148253248869,-1.8623430817467594,
0.0005218257781743,0.0046235194491761,0.0046683947419736,-0.0001574046868744,-0.0001382279541561,
-0.5017085745215428,-6.0467566606509120,3.6769898041979263,-1.0673454373145170,1.1362820763678718,
    ]);
    return a;
}

pub type AMat5 = na::SMatrix<f64, 5, 18>;

pub fn init_AMat5_pos() -> AMat5 {
    let a = AMat5::from_row_slice(&[
15.8312199273174397,-0.0120425882384158,-0.0746181509977484,48.4290700003549972,-0.4436603887830321,-0.6734649513605961,0.0087605247509007,0.0212400970627802,0.5067072548360047,0.4633429374413825,0.1027441582754757,-0.0044315680229212,-0.0009268311965047,0.0000436370323316,-0.0007339286727709,0.0000380322625124,-0.0203428316509982,-0.0042808555267044,
0.0080071798600724,8.2325235562609311,-0.0000308570315855,-0.2641902458930971,9.2392927505296978,0.0024761138796813,0.2545248700935174,-0.0470564000465378,-1.1177020670269124,-1.2264090198575879,-0.8429071656617971,-0.0442273234683109,0.0004114239998660,0.0002924202091842,-0.0000032769637049,0.0000092953807783,-0.1683707183940906,-0.0515500660810827,
0.0078434871290310,8.2334808580340919,-0.0000300902757251,-0.2644206228185210,9.2399294017111906,0.0024763042027292,0.1204155018236461,0.0225269819483873,-1.1698920305628109,-0.7441350055029123,-0.6741045711533067,-0.0097861541347817,0.0015334454080378,0.0000655427237964,0.0017953189112880,0.0002351741111090,-0.1278031738286389,0.0313962034337913,
0.0562729995457796,0.0004799324395161,3.4709451814268673,-0.3118082874550331,0.0031212112292790,3.0402984610150474,-0.1400046049052178,-0.0104202480691859,-0.1881564037443438,0.5245218395588421,0.1896088565955703,0.0065339712299302,-0.1341976139175835,-0.0222138505760713,-0.2131633439715113,-0.0254100971038089,0.0293566019807386,-0.0091078142224356,
0.0564488786584021,-0.0004181679101535,3.4709468540759363,-0.3115713722370168,0.0025184364204553,3.0403026546686118,0.1428356755743354,0.0101392571907842,0.2029600676603447,-0.5343610212159532,-0.2006655869388621,-0.0072949942362622,-0.1342046144971250,-0.0220417348418835,-0.2128603428714396,-0.0250756786511723,-0.0250347332722895,0.0096913059927722,
    ]);
    return a;
}
