

% =========================================================================
%
%                  IMU_kalman_filter
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 3月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.寻找高斯分布的特点
%--------------------------------------------------------------------------

clear;
close all;


l1 = randn(1,1e6) * (sqrt(2.5*0.25));
l2 = randn(1,1e6) * 3.5;
l3 = [12.6954000000000;12.6954000000000;12.6954000000000;12.6954000000000;12.6954000000000;12.6954000000000;13.5502000000000;13.5502000000000;12.8185000000000;3.04174000000000;3.03946000000000;3.01328000000000;3.06507000000000;3.07056000000000;24.0154000000000;24.0391000000000;13.5546000000000;12.6040000000000;5.01736000000000;7.48773000000000;12.7610000000000;6.28625000000000;12.8263000000000;13.9360000000000;16.3874000000000;15.1986000000000;6.34941000000000;6.30363000000000;13.8639000000000;10.1070000000000;10.1816000000000;11.4258000000000;13.8200000000000;-1.25423000000000;16.3035000000000;17.6916000000000;7.61953000000000;15.0873000000000;11.2862000000000;3.80072000000000;3.77250000000000;10.1600000000000;5.04230000000000;7.50559000000000;10.2072000000000;8.90819000000000;12.6368000000000;8.82162000000000;6.34984000000000;7.52422000000000;10.1106000000000;6.36404000000000;1.25085000000000;5.10480000000000;7.60080000000000;14.9699000000000;12.7937000000000;8.80123000000000;12.5436000000000;11.4579000000000;8.84670000000000;13.9667000000000;10.1087000000000;11.2557000000000;11.2880000000000;18.9095000000000;12.7949000000000;7.56278000000000;20.2202000000000;15.2372000000000;11.4300000000000;-2.27374000000000e-13;20.2668000000000;12.7624000000000;7.58175000000000;11.5099000000000;7.60114000000000;7.56360000000000;10.0301000000000;6.31738000000000;8.78067000000000;6.25556000000000;17.8654000000000;15.1259000000000;12.4785000000000;11.3409000000000;13.8030000000000;11.4280000000000;11.4561000000000;3.78138000000000;20.5134000000000;12.6970000000000;12.5734000000000;12.5708000000000;12.7606000000000;11.3448000000000;14.1096000000000;-1.26659000000000;6.26702000000000;6.31890000000000;6.34998000000000;6.27071000000000;13.9346000000000;8.93656000000000;8.82481000000000;17.7366000000000;12.7604000000000;10.2366000000000;12.5431000000000;17.9078000000000;11.4567000000000;10.1595000000000;11.4564000000000;10.1090000000000;13.7946000000000;13.9660000000000;16.5091000000000;15.1634000000000;10.1598000000000;2.53982000000000;11.3712000000000;5.11766000000000;11.3464000000000;15.2356000000000;10.1842000000000;8.91360000000000;13.8958000000000;8.86704000000000;10.1085000000000;7.59926000000000;8.82453000000000;8.88889000000000;11.5138000000000;12.6999000000000;6.34822000000000;13.9312000000000;8.82214000000000;8.75711000000000;17.5976000000000;15.3500000000000;12.6643000000000;5.04034000000000;8.93271000000000;11.2850000000000;7.63396000000000;8.81888000000000;11.4003000000000;5.06577000000000;20.2150000000000;11.5096000000000;11.3440000000000;13.8652000000000;10.1563000000000;11.3710000000000;10.0581000000000;11.4274000000000;12.4754000000000;8.79865000000000;11.1987000000000;7.48665000000000;2.55178000000000;5.02953000000000;10.0324000000000;12.5393000000000;12.6664000000000;6.34954000000000;12.7277000000000;10.2354000000000;1.25690000000000;12.7914000000000;13.8310000000000;11.2845000000000;5.02762000000000;7.65449000000000;11.3990000000000;10.0838000000000;11.2566000000000;7.54305000000000;6.39506000000000;12.7256000000000;11.3144000000000;7.52386000000000;8.84649000000000;11.4548000000000;20.1105000000000;6.37960000000000;8.82232000000000;8.84429000000000;6.25530000000000;15.1588000000000;10.0581000000000;13.8621000000000;12.6282000000000;8.90934000000000;10.0810000000000;11.4000000000000;8.95468000000000;12.6988000000000;12.5975000000000;13.7914000000000;13.7936000000000;8.82280000000000;8.86404000000000;13.7566000000000;5.06611000000000;11.3688000000000;16.5073000000000;8.84443000000000;8.82128000000000;18.7614000000000;8.80017000000000;12.6345000000000;6.31863000000000;5.04091000000000;10.1810000000000;2.52726000000000;6.28553000000000;7.65595000000000;1.26688000000000;1.26661000000000;17.9049000000000;7.58110000000000;7.56005000000000;3.79960000000000;8.88970000000000;10.0550000000000;13.7602000000000;1.27315000000000;12.6985000000000;10.0842000000000;10.1540000000000;12.5707000000000;3.81939000000000;16.5463000000000;7.59981000000000;12.3851000000000;12.6030000000000;2.53993000000000;12.6975000000000;13.9691000000000;5.09111000000000;6.24028000000000;12.5366000000000;15.1552000000000;6.28851000000000;11.3709000000000;8.82340000000000;8.84274000000000;0;11.4013000000000;7.62023000000000;15.0841000000000;14.1074000000000;8.88776000000000;8.86665000000000;13.8616000000000;13.8241000000000;8.97688000000000;7.56241000000000;11.3661000000000;7.48649000000000;12.7291000000000;11.3944000000000;10.1089000000000;10.0324000000000;3.79051000000000;8.77728000000000;8.95487000000000;17.5870000000000;7.52456000000000;17.8627000000000;10.0572000000000;6.39527000000000;7.58119000000000;12.6989000000000;15.3488000000000;11.3716000000000;6.39447000000000;5.02892000000000;6.36367000000000;12.7280000000000;12.5105000000000;16.4694000000000;8.73604000000000;17.7312000000000;7.59998000000000;13.8989000000000;12.5094000000000;7.65662000000000;10.1567000000000;7.67530000000000;6.25621000000000;12.8239000000000;17.6021000000000;-1.13687000000000e-13;7.65852000000000;8.86850000000000;7.56263000000000;11.5153000000000;8.88857000000000;10.1844000000000;12.6322000000000;15.1221000000000;11.3410000000000;8.84744000000000;21.7985000000000;6.33244000000000;11.4844000000000;15.2724000000000;11.5455000000000;8.77999000000000;8.77923000000000;11.4297000000000;11.2606000000000;7.56266000000000;13.8996000000000;12.6048000000000;3.76254000000000;7.65493000000000;13.8643000000000;11.3980000000000;15.0846000000000;7.48695000000000;8.86611000000000;16.4603000000000;15.1999000000000;15.2764000000000;13.8276000000000;8.86649000000000;13.8658000000000;16.6661000000000;7.61768000000000;16.4645000000000;12.5426000000000;14.0023000000000;3.77172000000000;3.78065000000000;12.5347000000000;10.0563000000000;7.61805000000000;10.1812000000000;15.1964000000000;12.6019000000000;19.0003000000000;13.9685000000000;14.0342000000000;10.0082000000000;11.3114000000000;3.79012000000000;11.3706000000000;6.33283000000000;7.54268000000000;16.6283000000000;10.2072000000000;10.1066000000000;13.9287000000000;16.4634000000000;10.1855000000000;11.3382000000000;14.0680000000000;11.4560000000000;13.9345000000000;7.59716000000000;11.5108000000000;6.26957000000000;11.2868000000000;12.6643000000000;3.78981000000000;17.8661000000000;10.1807000000000;15.1199000000000;7.56191000000000;15.0085000000000;5.00330000000000;8.95466000000000;10.0580000000000;8.86328000000000;10.1310000000000;7.48672000000000;5.04131000000000;10.1090000000000;6.34813000000000;8.91125000000000;7.59696000000000;11.4251000000000;21.5873000000000;10.1318000000000;11.4577000000000;11.4541000000000;5.11743000000000;8.84214000000000;5.02893000000000;3.79945000000000;10.0798000000000;16.4248000000000;13.9659000000000;12.6047000000000;10.1090000000000;7.65566000000000;2.53256000000000;11.4830000000000;11.3977000000000;13.9650000000000;10.1819000000000;11.2869000000000;5.10409000000000;8.84519000000000;11.4245000000000;15.2376000000000;10.1833000000000;11.3136000000000;15.2378000000000;10.0831000000000;16.4704000000000;16.5887000000000;7.63696000000000;6.42667000000000;8.82501000000000;11.5641000000000;10.1576000000000;13.8636000000000;10.2082000000000;10.2317000000000;1.26011000000000;10.1319000000000;7.52314000000000;2.53908000000000;20.3681000000000;11.4845000000000;6.34975000000000;5.05248000000000;7.65553000000000;5.11644000000000;15.2755000000000;18.9496000000000;8.90866000000000;16.6325000000000;5.05314000000000;13.8987000000000;3.80000000000000;12.6642000000000;8.90892000000000;10.1858000000000;7.58032000000000;16.3447000000000;-1.26334000000000;11.3129000000000;11.4303000000000;7.63825000000000;12.6364000000000;8.86500000000000;6.31954000000000;11.3405000000000;15.1251000000000;15.2725000000000;5.10439000000000;10.0827000000000;15.1986000000000;2.50183000000000;10.1584000000000;6.25539000000000;7.58048000000000;5.01632000000000;10.1062000000000;13.9328000000000;12.6363000000000;13.8634000000000;13.7220000000000;3.81890000000000;8.79777000000000;21.4268000000000;8.86730000000000;17.5120000000000;13.9330000000000;7.56326000000000;13.9677000000000;7.61870000000000;1.26665000000000;16.5884000000000;10.1069000000000;10.1334000000000;6.27083000000000;6.34950000000000;10.1068000000000;24.1263000000000;2.56517000000000;13.9269000000000;13.9327000000000;6.34848000000000;10.1862000000000;12.7311000000000;11.3447000000000;11.3755000000000;16.5115000000000;3.79059000000000;12.5436000000000;10.1068000000000;3.83757000000000;16.3874000000000;12.6656000000000;10.2076000000000;11.2311000000000;15.1998000000000;17.5976000000000;12.6984000000000;11.2320000000000;5.10537000000000;13.9674000000000;11.3146000000000;11.3975000000000;6.28512000000000;7.54266000000000;16.4604000000000;7.58135000000000;10.0278000000000;15.1236000000000;5.04179000000000;8.82230000000000;14.0655000000000;11.3141000000000;8.84386000000000;3.80935000000000;18.9480000000000;7.52092000000000;13.9998000000000;15.1638000000000;11.4262000000000;8.82242000000000;6.31729000000000;5.05333000000000;7.57852000000000;16.4577000000000;8.86433000000000;7.50459000000000;12.6023000000000;12.6333000000000;3.79138000000000;7.54200000000000;6.33115000000000;8.93170000000000;5.10360000000000;5.00399000000000;6.22248000000000;8.88776000000000;9.98138000000000;11.3407000000000;7.56157000000000;16.5004000000000;8.95215000000000;16.5048000000000;7.59672000000000;8.79936000000000;8.86703000000000;6.28736000000000;16.1771000000000;13.8984000000000;10.0818000000000;15.0481000000000;21.5218000000000;12.7571000000000;6.34815000000000;11.2528000000000;10.0519000000000;12.6269000000000;2.52095000000000;5.05340000000000;10.0552000000000;3.79699000000000;10.2074000000000;8.80084000000000;10.1530000000000;10.0292000000000;6.28680000000000;5.07800000000000;15.1571000000000;16.5081000000000;11.3692000000000;10.1327000000000;14.0340000000000;15.1159000000000;17.5502000000000;10.1849000000000;8.86510000000000;5.01609000000000;12.6338000000000;7.69449000000000;5.04162000000000;7.63685000000000;10.0294000000000;6.25317000000000;13.8624000000000;15.1969000000000;12.7568000000000;17.8597000000000;3.79963000000000;6.26851000000000;8.86461000000000;12.6641000000000;8.97653000000000;6.30100000000000;10.2092000000000;8.88668000000000;15.2737000000000;2.51387000000000;6.33351000000000;10.1347000000000;11.3730000000000;7.61853000000000;12.7905000000000;11.3983000000000;17.5544000000000;8.84385000000000;8.84566000000000;9.98325000000000;11.3987000000000;12.7294000000000;8.84062000000000;6.33397000000000;5.09214000000000;8.75579000000000;8.95315000000000;10.0343000000000;12.6669000000000;13.9336000000000;7.67596000000000;1.26378000000000;14.0342000000000;8.88969000000000;11.3693000000000;7.48573000000000;3.77298000000000;7.65520000000000;11.4038000000000;5.05403000000000;11.4562000000000;13.9665000000000;7.56439000000000;7.58029000000000;5.01792000000000;13.8627000000000;11.3710000000000;11.3936000000000;6.34741000000000;14.0027000000000;2.52053000000000;10.1848000000000;15.1600000000000;8.88671000000000;17.6424000000000;24.0650000000000;8.84081000000000;10.0838000000000;8.80184000000000;10.0301000000000;6.33152000000000;6.30119000000000;15.1565000000000;16.6244000000000;5.07943000000000;14.0026000000000;7.61836000000000;6.37938000000000;12.7936000000000;6.31673000000000;6.38029000000000;8.84582000000000;16.4690000000000;11.3430000000000;10.2088000000000;8.75705000000000;11.2889000000000;11.4568000000000;11.3713000000000;3.80852000000000;6.25283000000000;6.28624000000000;15.0838000000000;3.80012000000000;10.2083000000000;11.3436000000000;11.3128000000000;7.58299000000000;10.0840000000000;7.65642000000000;10.1097000000000;11.3437000000000;3.79039000000000;8.78134000000000;11.3455000000000;10.1324000000000;11.3985000000000;11.2862000000000;16.6302000000000;11.4005000000000;8.82403000000000;11.4537000000000;12.6068000000000;6.28475000000000;15.0496000000000;12.7326000000000;8.77859000000000;16.5445000000000;6.30203000000000;17.7360000000000;13.8616000000000;7.63675000000000;8.88830000000000;10.0785000000000;12.6308000000000;12.6971000000000;2.53911000000000;10.1337000000000;8.86342000000000;10.0600000000000;6.28597000000000;11.5122000000000;2.53311000000000;10.0825000000000;6.30046000000000;7.61929000000000;15.2768000000000;8.88751000000000;10.0564000000000;10.1318000000000;10.1062000000000;7.50529000000000;15.1975000000000;10.1559000000000;11.3407000000000;11.4823000000000;11.4850000000000;14.0017000000000;7.60016000000000;10.0561000000000;8.84418000000000;8.95158000000000;3.78979000000000;13.7905000000000;8.75657000000000;3.81780000000000;-2.27374000000000e-13;7.56249000000000;11.4535000000000;5.06529000000000;6.39528000000000;11.3168000000000;6.33273000000000;5.06734000000000;15.1539000000000;10.1845000000000;11.4273000000000;7.58177000000000;10.0826000000000;16.5476000000000;20.0568000000000;8.88550000000000;5.07869000000000;20.2574000000000;1.27907000000000;17.7738000000000;2.52123000000000;16.4169000000000;3.80952000000000;8.82168000000000;8.86512000000000;7.67387000000000;6.37976000000000;10.0336000000000;15.2706000000000;5.07854000000000;10.1838000000000;5.07862000000000;12.7570000000000;6.30109000000000;10.0772000000000;-6.82121000000000e-13;16.3817000000000;8.93220000000000;12.5682000000000;11.4255000000000;10.1347000000000;6.33338000000000;15.1550000000000;5.09112000000000;8.84414000000000;11.3656000000000;7.56189000000000;2.55807000000000;14.0676000000000;7.54205000000000;14.0025000000000];
zeta = 2;
count_perset(l1 ,zeta , 0 )
ksdensity(l1)

function out = count_perset(l ,zeta ,miu)%%数列、标准差、均值
cunt = 0;
for i = 1:length(l)
    if l(i)>miu-zeta&&l(i)<miu+zeta
        cunt = cunt+1;
    end
end
out = cunt/length(l);
end




