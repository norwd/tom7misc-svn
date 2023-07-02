// Generated file! Do not edit.
#line 1 "guitar-head.cc"

// Note: In guitar.cc, this file is assembled from guitar-head.cc,
// guitar-tail.cc, and some data generated by gencc.cc. Don't edit
// guitar.cc directly.

#include "guitar.h"

#include <string>
#include <vector>
#include <string_view>
#include <array>
#include <optional>
#include <unordered_map>

#include "util.h"
#include "base/logging.h"

using namespace std;

#line 1 "(generated)"
[[maybe_unused]] static constexpr int RADIX = 66; [[maybe_unused]] static int CharNum(char c) { if (c >= ',' && c <= '9') return c - ','; if (c >= 'A' && c <= 'Z') return (c - 'A') + 14; if (c >= 'a' && c <= 'z') return (c - 'a') + 14 + 26; return -1; }
static constexpr char DATA[] = "8,g,0,0/-.-,02220,,2225577655-0,0.-.00022105322,,577555.0,0.,./,0121,53,54,,,78,8/0,,././,01/1/5,454,,,787800,0--.0,0--000022005,-45510,00-..00223055,-355777552/,00-00,0--..,0020030,000..002030,,223357575540,0/22/,0122-,,7645,,78A950,0/..,,0322,54322,,,76656/57,,,,,02,,,577,,,70,0//.-,022225,465,5,767580,0//00,0--22544455,7767790,0/0.-,02020,,2223575655A0,,/0./,0102,5,564-,,7889c/00/0,,0./-.-002020B0,0/0,1,030215,566-,,7A89C1-0/-0-00/000545455575457,,7687D0,0/00/,0100-5454455656,7E0,0/001,0-02132322354546,F0,0/0/05453,,5,5656,,7686G0,0/01,,0201-575658,,7688H0,0/0..,0002054--33555655I0,0/00/,01020545445,,7889J0,0/022000022544453575675K000/---002120,,2224,,7999L0,0/1-/,0112,546644,,7899M0,0/.--,0312-5432--5876--N0,0----00/10,2222245,6657O0,0/--.,00---546-3-556655P0,0//-.,00122544454556675Q0,0./.0,02,125,4555577575S0,0.01,002010,,2213575555T0,0101,,,1213567585,,7888U0,0.000,0-010,35455575557R0,0./00,0--125344,,,77577V0,0.00.,00010535433555557W0,0.--,002110576555,,7998X0,0,11/,0111,5665,5,,7898Y0,0.10,,0-11053,4-5576557Z0,0.--.000110556557,,7798a0,0/-0-,0--0-54--5-,,7657b0,0.-00,0-21053-455,,7557d/-0/-.--0/-.0,42255e/,,0-.-.0/-.-,,0220g/00/-.-00222000/22,-c,0,10./.113331653336688766-0,1/./,1133216433,6688666.0,1/,/0,1232,64,65,,,89,9/0,,/0/0,120206,565,,,89890011331163356,688,68888AB810,11./,1133416,33466888662.,1131166886830,111//113141,,,34468686640,10-/.,12-30,1233,,,875650,111//11314165433,6,87766/68,,,,,13,,,688,,,70,100/,,133336553336,878680,1.0/.,10011655566,,878890,101/,,13131653334686766A0,10--.,1213,6,675,,,9AABB0,10/-.,14132654,-46,677,C01101116565,,686768,,8798D/,1011065655667,7-8E0,101126545-,65657,696768F0,10101,1-1346564-,,,6767G0,101/-,1012,,56666686769H0,10--1,12131656556676766I0,0/--0,12131656556676766J0,100-/111133655564666788K0,10...113231,,,3356,776,L0,10230,1223,657755,,89AAM0.10/..,1423,,,8775,,8BAAN0,1....,1021,657565,,87A8O0,102//11123165,,45667766P0,1.0..,11233655565667786Q0,1/0/1,,33236,566,688686S0,13121,,3324686666,,8A99T0,1212,,,23246,665,,,8999U0,1/111,33324,46566686668R0,1.0/-,1/01,6455,-,88688V0,1/1//646544666668,,8899W0,1/..,113221687666,88AA9X0,1/-.-,1222,6776,6,,89A9Y0,1....,1/21-11321-687668Z0,1/21/111221667668,,88A9a0,10.1.,10,1165,56,,,8768b0,1/.1,,1331-64,56-,,8668.h,0,,-/0/,21/0/,24442799877-0,,-/0.224432,5443,799777.0,,-.,.,20,01,2343,75,76,/0,,-.-.,23131,5-76478979700,,-/0-,,//02224422,4447710,,-/00,2--022244527999772/,22422,2//0077997930,,-/.0,2220,22425279797740,,-.0/,21,01,2344,,8-86750,,-00/,2100,76544,7,988,6/79,,,,,24,,,799,,,70,,-/-/,2110,,2444446447480,21/--,21122766677,9989990,,-/./,2120,224242797877A0,,-../,,1201,2324,7,786,B0,,-0./,212,3,252437,788,C0221222,4-44576767,797879D0,21221,2324-7676677878,9E0,21223,4-24376768,7,788-F0,,-21/,21212,3-2427675,,G0,2123,---745,5444579787AH0,,--./,21200,22242,64455I0,,-../,21221,6-465767667J0,,-1./22224476-6-57,789,K0,,.///,21///224342,,4446L0,,-.//,2334,7688667888,7M0,,-0///210//,2-343,,-886N0,2////,2131,,6-6467,8879O0,,--//222342,6-456778877P0,,-1//,211/-,22344766676Q0,,-/-.,20102,24,347,6777S0,,-/..,24232,,4435797777T0,,-...,20201,2323,,5-765U0.--/.-,20222,2423-797779R0,20/--,2012-,4-434,99799V0,,....,20220757655777779W0,,-//.,20//-224332798777X0,,-./.12333,,2333,7887,7Y0,2032-,4-332,5-676798779Z0,,--/.,20320222442778779a0,21/2/,2444-,,-647,,9879b0,20/0-,20/2,,2443-,,9779f//--/0//21/0/,,1/02i/,--/0/22444222144,j/,.-/0/,.1/0/32444,k/,/-/0/,/1/0/,/1/02l/,0-/0/,01/0/52444,/c,0,,.010,32010,35553,,5558-0,,.01/,,101/3355438AA888.0,,./,/,31,12,3454,86,87,/0,,././,34242,,45458,787,00...01.33553385578,8AA,8A10,,.011335563,,556,8AAA882/,33533,3001188AA8A30,,.0/1,333113353638A8A8840,,./10,32-12,34-5,87,97850,,2110,3211,87655,8,A99,6/8A,,,,,35,,,8AA,,,70,,.0.0,3221,,355558,798,80,,.-..,00010,3223387778890,,.0/0,35353,,55568A8988A0,,.//0,,2312,343538,897,B0,,.1/0,36354878--68,899,C0,,.-/.,323338787,,8A898AD0,32332,34-368787788989,AE0010110,3233487,7-687879,F0,,.-/-,32323,3535-8786,,G0,,.-//,3234,,655568788,,H0..../0,32311333353888988I0,,.//0,32332,3435,878778J0,32355333355877-668889AAK0,..000,32000335453,,5557L0,,./00,344538799778999,8M0032100,324,4,3645,8,9998N0,30000,3243,,55557,,998AO0,,..00,32411333453889988P0,0.-0,,32200,33455877787Q0,..0./,31213,,55458AA8A8S0,,.0//335343,,55468A8888T0,,.///,3434,,,45468,887,U0,31333,,5746,687888A888AR0/,.0..,3123,867766,,A8AAV0,,..//,3133186876688888AW0,,.00/,3100,,354438A9888X0,,./0/23444,,3444,8998,8Y0,3143,3654338697,,8A988AZ0,...0/,3143133344388988Aa0,32030,32-3387,78,,,A98Ab0,,101.,3103,,31,33,,A88A0s,0-//.--,,/121,43121,46664-0-//----//120,466549766,,.0,,/0,0,42,23,4565,97,98,/0-./-/-,,/0/0,4535,9,898,00///12/446644-666--9BB,9B10-///---//122446674,,667-20--//-/-///-/,44644,4112230-/-/--,44422446474,,667740-./.,,,,/021-43,23-4566-50-0/..-,4322,,4766,98766-609B,,,,,46,,,9BB,,,-//,,,70-//./-.//1/1,4332,,4666680-//.//,43344,66666-8889990-/-.--,4342,446464,,6667A0-.-.0-,,/001,,3423,4546,B0-0-..-,,/201,474659,9AA,C0-/-.-/1,/10/443444-66667D0-./.0/-21001,43443-6546-E0-0-.0/-21221,43445-67667F0-/-.-.-2-131,434349897,,G0-/-.-0-2--01,4345,-76667H0---.----1101,43422,44464I0-.-.---2-021,43443989889J0-/-./----.//-2332-444466K0-/..--,,/111,43111446564L0-...1-,,/011,4556,-8AA88M0-0/.11143211,4756,-87668N0-/..-/1//11/,4354,,,6868O0--..----3112444564-86678P0-/..//,43311,44566988898Q0-//-/--//1/0,42324,,6656S1-//-0--/-----//100446454,,6657T0-./000,4545,,,56579,998,U0-/---/-//-0/,42444-79897R0-//-//0,/1//,42344-7889-V0-----/--210-,42442979877W0-/.----//110,4211,446554X0-..-,-,,/01034555,,4555,Y0-/.--/-,1110-7654-9BA99BZ0--.--/,///10,42542444554a0-//.-/,43141,43,44,,BA9Bb0,,0..0,5325,,53-35,57-65k/,//.--,//121446664n0,//---,//--0,//120446654m/,13121,13124686699p/,1/---,12124676999q/,,-.--,23121,2/121r/,,----,2212,,,9999s/,,..--,33121886699t/,,.---,3212,,3665,e/.//.--,,0.--,,0121f////.--,,1121///121g/0//.--0//121,,2121h0,,,.--1//.--1431211//12,v/,,0---.//---.//--0w////---///--0///120x/0//---0//--00//12-y/1//---14212,,,31201j,0.00/..,,0232,54232,57775-0.00...,,0231,57765A877,,.0,,01,1,53,34,5676,A8,A9,/0.,-.-.,,0101456464A,9A9,00.00,.0000230557755A779A,10.000..,,0233557785,,778,2/..00.0,55755,5223330.0.0..355533557585,,778840./0/-,,,01-2,5677,A9,B9A50,,0//.,5433,,5,776A9877,6/.0,,,,,57,,,.00,,,70.,0/0.,00202,5443,55777780.---..,,0/00,22232,5445590.0./..,00212557575,,7778A0.-./-.,,0112,,4534,5657,B0.-.//,,,0312,58576A9A,B,C0.0./.0,,0/10554555,77778D0.-.--.,54554,56-78A9A--AE0.-.-/.232332,54556A9A9B,F0.0././,,0/1/,54545A9A8,,G0.0./.1,,0/11,5456,,87778H0.../..000012,54533,55575I0.-.--.,,0112,54554,56575J0.0./0..../00,54577555577K0,,0/.-.0//..,00222557675L0.-//--,,0122,5667,A9BB99M0.-///-254322,5867,,,CBB9N0.-/-.-,,//.0,5465,,,7979O0..//..,,0022555675A9B-8-P0.---.-,,0/0-,54422,55677Q0.,-...,,0201,53435,,7767S0.0....,,0211557565,,7768T0.,..-,,,0111,5656,A8AA-,U0.0...0,00-11,53555A8A-AAR0.00.001,0200,5345,A89-A,V0.....0,,0011,53553A8A988W0.0/...,,0221557665,,CAA9X0.//.--,,012145666,,5666,Y0.0/..0,,0-21,5365,A8,9A-Z0../..0,00021,53653555665a0,,0/.0,,0-32,54252,54-55b0,,0..0,5325,,53-35,57-65l/,00/..,00232557775q.,24232,24235s/,,./..,34232,30232d/--0/..,,//..-00/..g/0-0/..000/,,000232i0,-0/..200/,,25423220023-o/,00...,00..15577652c,0/110//,11343365343668886-0/11///,11342,,4342668876.0/-,/.,,,12,2,64,45,6787,/0/,./.,/01/1/,,121256757500/,,.//111341668866B88AB,10/111//,11344668896,,889B2///11/1,66866,6334430/1/1//,,1324,6664466869640,,10.//010,/,,1243,6788650,,100/,,100,,6544,BA988,6//1,,,,,68,,,/11,,,70/,.0/,/,101/,11313,6888880/...//,,1011,33343,6556690/1/0//,11323,6564,668686A0/,/0.,,,1223,45645,6768-B0/,/00,,,1423,6544-,69687C0/1/0/1,,1021665666,88889D0/./..//./..-,01021,65665E0/./.0-,6546-,65667B,9A8-F0/./-/-,,1020,65656BAB9,,G0/./////1/0/2,,1022,6567,H0/./.--111123,65644,66686I0/./../,,1223,65665,67686J0//.0--///011,6554-666688K0/100//,11333,65333668786L0/.00..,,1233,65785,6778,M0/,000,365433,6978,,,988AN0/.0././,00/1,,103166576,O0//00//,,1133,657-6666786P0/.../.,31133,65533,66788Q0/,.///,11312,64546,6857,S0/1////,,1322668676,,8879T0/-//.-,,1222,6767,,,7879U0/-/./-/1///1,64666,,8A79R0/-..///11/11211311,6456,V0/-/.--/////1,,1122,64664W0/10///,11332,6433,668776X0/00/,/,,123256777,,6777,Y0/-0././-0./,/10//1,6476,Z0//0//1,11132,64764666776a0/.,.//,,10/1,65363,65,66b0,,1//1/11//1,6436,,64,663h,00/---0022100,,2454476454-00.--00022000,,2453779987.00.,0/,,,23,3,75,56,7898,/00.,0/-0,/0/-012020,,2323000---00222452,744577799771000--.0022200,224557799A72100-/.,00-/0,0--/.0002202,779773000....020200,777557797A7400/,-/0,,23-4,76--6,7899,500/.--,0,211,,,211,,7655,6/02,,,,,79,,,022,,,700/----0,/10,,22424,79999800-----0///00,22122,76677900/---.020100,22434779797A00,01/,,,2334,767-6,7879,B00/.--.0,011,,,2534,7A798C00----.0/0/0,020102776777D00/0//00101-2,76776,98--AE00-.--.0/0/1,,76778,9A-9AF00/-.0.0,0101,,2131,76767G00/-0-.020133,22133,7678,H00/--..222234,76755,77797I00/0//0,,2334,76776,78797J00-/--.0,0122020120,76779K00/---/021100,22444,79897L00/11//0111,,,,2344,7889,M0,,.--/0,111,476544,7A89,N00----/0/1/0/,,1102,7687,O00/--./001100,,6-45777897P00///0/001120,76644,77899Q00,/000022020,22423,75657S0020000,22433,75757779787T00,,0/.0,00/,,,2333,7878,U00--00.020002,,4-33,75777R00.-/0-,22022322422,75677V00,00.,000002,,2244,75775W00.--0/021000,22443,79887X00110,0,,234367888,,7888,Y00--00/021002,,4-43,7587,Z0001002,22243,75875777887a00--/-0,,0/.1,76474,76-77b00.-/00,,2002,,4-53,7547,f///---0//--00,,1100k/,/--00,/---0,/2100q/,,---0,22100,22454d/-/--00-/---0,,/100e/./--00./---0,,01004b,010...,133211,33565587565-0133111,,3564,,656488AA98.01/,10,,,34,4,86,67,89A9,/0,,-.-.1,0101,,3434789797001,,011133,1333356388AA8810,,../1133311,3356688AAB830,,..//131311,8886688A8B840,,3201,,3465,87AA7,89AA,5010/..,1,322,,,322,,8766,6/13,,,,,8A,,,133,,,70,0....1,021,,33535,8776,80,.....100011,,3233,8778890,,.../131211,3354588A8A8A01,120,,,3445,67867,898A,B01,122-,,3645,8786-,8B8A9C010101,131213,,3243887888D010-0,/1010011212,3,87887E0/./../10102-,,1223,8788-F0,-.../,,1212,,3242,87878G0101111131244,,3244,8789,H0111211333345,87866,888A8I01,120,,,3445,87887,898A8J01.0../131231111233,878AAK0132211,33555,8796,88A9A8L01022001222,,333455,899A,M010/-.-102-2-587655,8B9A,N0,....010201013221388798,O010.-/,112211,,33558889A8P0100010112231,87755,889AAQ01,011,133131,33534,86768S0131111,33544,8686,88A898T0,,-.-/1,110,,,3444,8989,U01...-/131113465564,86888R01,0113133133433533,86788V01/10//111113,,3344,86886W0,/..10132111,33554,8A998X01221,1,,345478999,,8999,Y01,20-1132113333554,8698,Z0112113,33354,86986888998a010,01,,,3213,87585,87,88b01/,01,,,3113133113,8658,5o,0,-///-,////2244322,-4676-0,-//.-,//222244222,-4675.0,-./.,20,10,,,45,5,97,78/0,-././2,121,234242,-454500,-//--,-/1--444674,-66--10,-//0-,--,0-244422,446772/,--/--22442422/1,,30,-/-0-242422,-4657,9997740,-.//,,-4312,-8576,9ABB,50,-0//.210//,2,4332,9877,6/24,,,,,-/,,,244,,,70,-////2,132,24,342,9887,80,-11//211122,44344,9889990,-/-/-,.///0242322,-4656A0,-.//0,-2312,,4556,78978B0,-0-/.,-0//0,22332,,4756C021/---,-/1/0242324,98999D0,-.1/0212112,-6556,98998E0,-01/0210--2,-2334,989-AF0,-/0/0242323,,4353,98989G0242355,,4355,-7656,989A,H0,---/-21---0222322444456I0,-.-/-212112,,4556,98998J0,-/-//211-0-242342,989BBK0,-/./-,.///1243322,-4666L0,-../1,-.//12333,,,-4566M0,-0./.,-0//1,-3332698766N0,-/1/121312-,-3324998A9,O0,--./-,--//1223324766666P0,-/.//211121,-3344,98866Q0,-//./2,1222244242,44645S0,-/-.-,-//.0242222,44655T0,-.-.,2,221,234252,,4555U0,-/1.02422245-6-5-,97999R0,-12--244244544644,97899V0,---.-202100222224,97997W0,-/..-243222,44665,9766,X0,-...12332,2,,456589AAA,Y02031--,-32-1243224,-6665Z0,--..-223224,-4465,97A97a0,-/1/-,,4324,98696,98,99b0,-/1.-,,4224,,46-5,9769,m/,1//22,1///2,14322d/--///-,44322,,///2e/.-///-,,0//-.1///,f//-///-/1//22/1///2g/0-///-,,232201//22h/11//2211///2,,3322o/,0//.-,0/222,0//22u/--//.--0//.-,44222v/.-//.-,,0/.-.0//.,w//-//.-,,1222/0/222x/,,22220-//.-00/222y.10/22210//226c,0,.000.32000,355433,55787-0,.00/.355333,,5333,,5786.0,./0/,31,32,,,56,6,A8,89/0,./-/-,,/0/0345353,,565600..00..30023,555785,A778A10,.001.,,0013355533,557882/,..0..3355353302,,30,.0.1.,,0011353533,,576840,./00-34-42-34-4,3,,568750,.100,32100,,,5443,,544,6/35,,,,,.0,,,355,,,70,.000032--3335,453,A887,80,.--..322233,55455,A99AA90,.0.0.353433,55767,A9A8,A0,./.0,3,342-,,5667,A9A8-B0,.1.0/3,344,,,5867,A9A,BC0,.-...32323,353435,A9AAAD0,.-..-3,-221323223,A9AA-E0,.-../10100132324,3,3445F0,.-.-.3231,,3,3434,A9A9AG0,.-./,353466,,5466,A9AB,H0,...0.321211555567,A9A88I0,.-..-32322,,,5667,A9AA9J0,.-.0032--11333455,A9ACCK0,.0/0.,,0002354433,55777L0,.//0,3244223444,-,,5677M0,.-/0/,,10023,444,7A9877N0..-/.,30000232423,354435O0,../0.3,-212334433,,5577P0,../00322232334453,A9977Q0,.0,/0,1203,3,2333355353S0,.0./.,,00/1353333,55766T0,././,,,/0/13,332,,,5666U0,,02/131,231353335,A8AA,R03,23353,53556,5755,A89A,V0313211333335,,5566,A8AA8W0,.0//.354333,55776,A877,X0,.///-34,3223443,3,,5676Y031,2323142,,354335555776Z0,..//.334335,55576,A8BA8a0,.-0..,,5435,,54,5,A97A7b0,100.,31,23,,,5335,A87A,7c,0//111/,,1114466544,66898-0//110/466444,,6897,,9897.0,/-,-.,/010,42,43,,,67,7/0,/0.0.,,01014,343,45646400//11//4,,344666896,B889B10//112/,,1124466644,668992/,//1//4466464413,,30,///--//1/2/,,112246464440,/011,,,65344565-,,,67-850,/.--,,,21104,655,,BA9996/46,,,,,/1,,,466,,,70,/..-,,/11114,656,,,686880,/..//433344,66566,B8A9890,/./-///1/1/,,1112464544A0,/./-.,/0/1,4,453,,,6778B0,/./-0,/2/104,455,,,6978C0,/.///411312464546,,6576D0,/.//.4343344545,6,BABBAE0,/.//021211243435,898998F0,/././4342,,4,4545,,6575G0,/./0,43-4-2464577,,6577H0,/./--////1/444544666678I0,/.//.,/0/1/4,453,,,6778J0,/./11////11413112444566K0//101/,,1113465544,66888L0,/001,4355334555-,,,6788M0,/.0-0,/201,435-5,8BA988N0//.0/,,,131343534,,,6586O0,/.0--///01/445546,,6688P0,//011433343445564,BAA88Q0//-.-/,,11014,344,466464S0//1/0/,,1102464444,66877T0,/0/0,4,443,456474,,6777U0,/-///,11102464446798897R0,/-.//,233444664667,6866V0,/-//-424322444446,,6677W0,/-0-///100/465444,66887X0,/000,45,4334554,4,,6787Y0,/-0//4253-,465446,B9CB,Z0,/-0/-///00/445446,66687a0,/.,//43,344,,6546,BA8B8b0,211/,42-34,,,6446,B98B,";
#line 1 "guitar-tail.cc"

// Note: In guitar.cc, this file is assembled from guitar-head.cc,
// guitar-tail.cc, and some data generated by gencc.cc. Don't edit
// guitar.cc directly.

#undef DEBUG_GUITAR

// Normalize the guitar chord name if possible.
// Sharps/flats on the base chord are rewritten to one of
// C#, Eb, F#, Ab, Bb.
static string NormalizeBase(string s) {
  if (Util::StartsWith(s, "Db")) {
    s[0] = 'C';
    s[1] = '#';
  } else if (Util::StartsWith(s, "D#")) {
    s[0] = 'E';
    s[1] = 'b';
  } else if (Util::StartsWith(s, "Gb")) {
    s[0] = 'F';
    s[1] = '#';
  } else if (Util::StartsWith(s, "G#")) {
    s[0] = 'A';
    s[1] = 'b';
  } else if (Util::StartsWith(s, "A#")) {
    s[0] = 'B';
    s[1] = 'b';
  }
  return s;
}

static inline Guitar::Chord ChordOfUnchecked(int base_num,
                                             int suffix_num) {
  return (base_num << 8) | suffix_num;
}

static inline pair<int, int> UnChord(Guitar::Chord c) {
  const int base = (c >> 8) & 255;
  const int suf = c & 255;
  return make_pair(base, suf);
}

Guitar::Chord Guitar::ChordOf(int b, int s) {
  CHECK(b >= 0 && b < NUM_BASES) << b;
  CHECK(s >= 0 && s < NUM_SUFFIXES) << s;
  return ChordOfUnchecked(b, s);
}

static string NormalizeSuffix(const string &s) {
  if (s == "major" || s == "maj") return "";
  if (s == "minor" || s == "min") return "m";
  if (s == "6/9" || s == "6add9") return "69";
  if (s == "m6/9" || s == "m6add9") return "m69";  
  return s;
}

int Guitar::BaseNum(string_view base) {
  for (int i = 0; i < NUM_BASES; i++) {
    if (BASES[i] == base) return i;
  }
  return -1;
}

int Guitar::SuffixNum(string_view suffix) {
  for (int i = 0; i < NUM_SUFFIXES; i++) {
    if (SUFFIXES[i] == suffix) return i;
  }
  return -1;
}

std::string Guitar::ChordString(Chord c) {
  const auto [base, suf] = UnChord(c);
  CHECK(base >= 0 && base < Guitar::NUM_BASES) << "Invalid chord " << c;
  CHECK(suf >= 0 && suf < Guitar::NUM_SUFFIXES) << "Invalid chord " << suf;
  return (string)BASES[base] + (string)SUFFIXES[suf];
}

std::string Guitar::FingeringString(Fingering fing) {
  const auto [a, b, c, d, e, f] = fing;
  auto Char = [](int x) -> char {
      if (x < 0) return 'x';
      if (x <= 9) return '0' + x;
      return 'a' + (x - 10);
    };
  char ret[7];
  ret[0] = Char(a);
  ret[1] = Char(b);
  ret[2] = Char(c);
  ret[3] = Char(d);
  ret[4] = Char(e);
  ret[5] = Char(f);
  ret[6] = 0;
  return &ret[0];
}

std::optional<Guitar::Chord> Guitar::Parse(std::string_view sv) {
  string s = NormalizeBase((string)sv);

  // Longer matches have to come first here, since
  // we use a totally greedy strategy!
  for (string_view pfx : { "C#"sv, "C"sv, "D"sv, "Eb"sv,
        "E"sv, "F#"sv, "F"sv, "G"sv, "Ab"sv, "A"sv, "Bb"sv, "B"sv }) {
    if (Util::TryStripPrefix(pfx, &s)) {
      // Matched.
      const int base_num = BaseNum(pfx);
      s = NormalizeSuffix(s);
      const int suffix_num = SuffixNum(s);
      // Invalid. (Note: Shouldn't continue looping since we've modified s.)
      if (base_num < 0 || suffix_num < 0) return {};
      const Chord c = ChordOfUnchecked(base_num, suffix_num);
      return {c};
    }
  }
  
  // No match for base.
  return {};  
}

// int8_t packed into lowest 6 bytes.
using PackedFingering = uint64_t;
static constexpr Guitar::Fingering Unpack(PackedFingering pf) {
  return make_tuple((int)(int8_t)(255 & (pf >> (5 * 8))),
                    (int)(int8_t)(255 & (pf >> (4 * 8))),
                    (int)(int8_t)(255 & (pf >> (3 * 8))),
                    (int)(int8_t)(255 & (pf >> (2 * 8))),
                    (int)(int8_t)(255 & (pf >> (1 * 8))),
                    (int)(int8_t)(255 & (pf >> (0 * 8))));
}

static constexpr PackedFingering Pack(Guitar::Fingering f) {
  return
    ((255 & (uint64_t)(int8_t)(std::get<0>(f))) << (5 * 8)) |
    ((255 & (uint64_t)(int8_t)(std::get<1>(f))) << (4 * 8)) |
    ((255 & (uint64_t)(int8_t)(std::get<2>(f))) << (3 * 8)) |
    ((255 & (uint64_t)(int8_t)(std::get<3>(f))) << (2 * 8)) |
    ((255 & (uint64_t)(int8_t)(std::get<4>(f))) << (1 * 8)) |
    ((255 & (uint64_t)(int8_t)(std::get<5>(f))) << (0 * 8));
}

// Sanity check the above for -1 especially.
static_assert(Unpack(Pack(make_tuple(-1, 3, 2, -1, 1, 0))) ==
              make_tuple(-1, 3, 2, -1, 1, 0));

namespace {
// Singleton parsed database. Depends on DATA symbol generated by
// gencc.
struct DB {
  DB() {
    data.resize(Guitar::NUM_BASES);
    int idx = 0;
    auto GetNum = [&idx]() { return CharNum(DATA[idx++]); };
    const int bases_in_data = GetNum();
    CHECK_EQ(bases_in_data, Guitar::NUM_BASES) << "Expected to be dense.";

    [[maybe_unused]] int ambiguous = 0;
    for (int b = 0; b < Guitar::NUM_BASES; b++) {
      // Not guaranteed to be in the same order.
      // This will be the index into the vector.
      const int base = GetNum();
      std::unordered_map<int, vector<PackedFingering>> &base_row = data[base];
      
      const int nsuf = GetNum();
      for (int s = 0; s < nsuf; s++) {
        // Suffix code (key into map).
        const int suffix = GetNum();
        vector<PackedFingering> &fingerings = base_row[suffix];

        const Guitar::Chord chord = ChordOfUnchecked(base, suffix);
        
        // Number of fingerings.
        const int nfing = GetNum();
        for (int f = 0; f < nfing; f++) {
          // Each is encoded +1 so that -1 is nonnegative.
          const int f0 = GetNum() - 1;
          const int f1 = GetNum() - 1;
          const int f2 = GetNum() - 1;
          const int f3 = GetNum() - 1;
          const int f4 = GetNum() - 1;
          const int f5 = GetNum() - 1;
          const Guitar::Fingering fingering =
            make_tuple(f0, f1, f2, f3, f4, f5);
          const PackedFingering pf = Pack(fingering);
          fingerings.push_back(pf);
          // There are some ambiguous fingerings in the upstream data,
          // so we just keep the first one we saw.
          [[maybe_unused]] bool inserted = rev.emplace(pf, chord).second;
          #ifdef DEBUG_GUITAR
          if (!inserted) {
            ambiguous++;
            /*
            LOG(ERROR) <<
            "duplicate fingering: " <<
            Guitar::FingeringString(fingering) <<
            "\nfor: " << Guitar::ChordString(chord) <<
            "\nand existing: " << Guitar::ChordString(rev[pf]);
            */
            printf("%s for %s and %s\n",
                   Guitar::FingeringString(fingering).c_str(),
                   Guitar::ChordString(chord).c_str(),
                   Guitar::ChordString(rev[pf]).c_str());
          }
          #endif
        }
      }
    }
    #ifdef DEBUG_GUITAR
    printf("There are %d ambiguous fingerings!!\n", ambiguous);
    #endif
    
    // sizeof DATA includes nul terminating byte.
    CHECK_EQ(idx + 1, sizeof DATA) << (idx + 1) << " " << (sizeof DATA);
  }

  // Outer vector is base chords, dense.
  // Then suffix, then each known fingering.
  vector<std::unordered_map<int, vector<PackedFingering>>> data;
  // PERF: These data structures don't depend on each other, so
  // we could perhaps only load the one we need (for embedded
  // applications that only call GetFingerings xor NameFingering).
  std::unordered_map<PackedFingering, Guitar::Chord> rev;
};
}

static const DB *GetDB() {
  static const DB *db = new DB;
  return db;
}

std::vector<Guitar::Fingering> Guitar::GetFingerings(Chord c) {
  const DB &db = *GetDB();

  const auto [base, suf] = UnChord(c);
  CHECK(base >= 0 && base < Guitar::NUM_BASES) << "Invalid chord " << c;
  CHECK(suf >= 0 && suf < Guitar::NUM_SUFFIXES) << "Invalid chord " << suf;
  const auto &m = db.data[base];
  auto sit = m.find(suf);
  // No fingerings known.
  if (sit == m.end()) return {};
  // Convert to external representation.
  vector<Fingering> ret;
  ret.reserve(sit->second.size());
  for (const PackedFingering pf : sit->second)
    ret.push_back(Unpack(pf));
  return ret;
}

std::optional<Guitar::Chord> Guitar::NameFingering(Fingering f) {
  const DB &db = *GetDB();

  PackedFingering packed = Pack(f);
  auto it = db.rev.find(packed);
  if (it == db.rev.end()) return {};
  return {it->second};
}
