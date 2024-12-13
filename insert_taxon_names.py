#Script by August Søndergaard

import re

# Path to the map file
map_file_path = "C:\\Users\\augus\\SKOLe\\5. semester\\Bachelor\\names.map"

# Newick string input (replace with your actual string)
newick_string = """(((Dypsis_lantzeana_1:0.002827,(0034:0.001470,Dypsis_lantzeana:0.001535)0.999021:0.001160)0.996562:0.001042,(((0037:0.003095,(0172:0.003908,((0179:0.003167,0069:0.002976)0.505006:0.000740,(0158:0.005373,(((((0036:0.002059,((0070:0.004183,0169:0.005362)0.642779:0.000009,0146:0.001180)0.884956:0.000187)0.999978:0.000621,0083:0.002556)0.372685:0.000623,(0033:0.001685,((0068:0.001657,0163:0.001212)0.867906:0.000291,0046:0.002095)0.605987:0.000127)0.988357:0.001655)0.757158:0.000722,((0007:0.003563,0005:0.002593)0.385824:0.000051,(0184:0.002236,0099:0.002397)0.912954:0.000625)0.460466:0.000420)0.784549:0.000213,(((0020:0.001553,0097:0.001403)0.839371:0.001449,0098:0.003446)0.544268:0.000067,(0048:0.003381,(((0090:0.002486,0013:0.002418)0.964129:0.000172,0029:0.003284)0.601931:0.000069,0073:0.002537)0.765152:0.000434)0.762062:0.000115)0.762163:0.000936)0.481899:0.000891)0.615585:0.000135)0.427340:0.000070)0.639354:0.000090)0.564300:0.000922,0050:0.003551)0.551098:0.000082,((((((((((((((((((((((2015:0.000487,((0102:0.000607,0103:0.000602)1.000000:0.001649,(0171:0.001152,(0174:0.001754,0095:0.000615)0.763458:0.001635)0.698948:0.000112)0.830874:0.000180)0.967105:0.002307,0053:0.002616)0.999429:0.000309,((0003:0.001954,0059:0.001741)0.948712:0.001392,(0120:0.000736,2014:0.000275)1.000000:0.002243)0.874771:0.001332)0.389932:0.000055,0123:0.002794)0.688963:0.000082,(((0067:0.001738,(0052:0.001055,0008:0.002049)0.972723:0.000654)1.000000:0.000940,(0072:0.001933,0183:0.001806)0.992453:0.000891)0.993905:0.000879,(0195:0.001784,0112:0.001017)1.000000:0.001175)0.398108:0.000054)0.975716:0.000953,(0111:0.001249,(0011:0.001579,0028:0.001678)0.992502:0.000282)0.999999:0.001495)1.000000:0.000771,(((0004:0.000870,0051:0.001136)0.999999:0.001499,((0202:0.000538,2012:0.001376)0.999999:0.001871,0006:0.002312)0.507910:0.000067)0.999996:0.000739,((0002:0.000941,0115:0.000623)1.000000:0.001100,0092:0.001932)0.999644:0.000779)0.999756:0.001347)0.837859:0.000190,(((0061:0.002303,0010:0.003187)0.999714:0.000343,0044:0.002977)1.000000:0.001779,0082:0.004932)0.999889:0.000135)1.000000:0.001792,(((((0149:0.000066,2034:0.000013)0.543656:0.003774,2052:0.002411)0.737676:0.001066,0045:0.002796)0.469734:0.003465,2051:0.002646)0.837868:0.003777,((0113:0.004670,((2033:0.000656,0204:0.001228)1.000000:0.005051,0091:0.002880)0.999596:0.000562)0.999999:0.001190,((0063:0.001427,0057:0.001275)1.000000:0.001278,0071:0.003005)1.000000:0.001257)0.769702:0.000296)1.000000:0.001705)0.943946:0.000621,((0012:0.000844,0021:0.000852)1.000000:0.002564,0181:0.004931)0.775400:0.000457)1.000000:0.001291,(((0109:0.002376,(0040:0.002147,0105:0.001954)1.000000:0.000593)1.000000:0.000908,0041:0.004585)1.000000:0.003102,(0078:0.004793,(0014:0.001607,0150:0.002024)1.000000:0.003144)0.999999:0.001059)1.000000:0.001114)1.000000:0.003294,(((((0055:0.001190,(0170:0.001375,0110:0.000427)1.000000:0.002890)1.000000:0.002416,(((0026:0.001069,0128:0.000599)1.000000:0.001061,0088:0.002484)1.000000:0.001123,((0074:0.001089,0116:0.001190)1.000000:0.000995,0009:0.002399)1.000000:0.001655)0.998166:0.000620)1.000000:0.005050,0031:0.010883)1.000000:0.001683,(((Loxococcus_rupicola:0.001904,1011:0.001554)1.000000:0.022657,(0187:0.001972,0077:0.002895)1.000000:0.008237)1.000000:0.001719,(0191:0.002604,0185:0.003891)1.000000:0.005019)0.712057:0.000324)1.000000:0.001518,((0196:0.002144,Marojejya_insignis:0.001051)1.000000:0.001843,0205:0.001820)1.000000:0.008019)0.653087:0.000773)1.000000:0.008655,0060:0.006904)0.760878:0.001270,((((((((0080:0.002345,0203:0.000180)0.695968:0.000051,0122:0.002391)0.780785:0.001109,0155:0.004066)0.922243:0.001032,(((0038:0.002191,(0124:0.000777,0081:0.000402)1.000000:0.002009)1.000000:0.001540,0160:0.003792)0.834646:0.000527,0027:0.002718)0.544128:0.000959)0.638983:0.001430,0087:0.005056)0.740993:0.000148,((((0106:0.001333,0126:0.000979)0.999999:0.001389,((0019:0.001487,(0153:0.003722,0201:0.002087)0.999900:0.000711)0.781963:0.000615,(0114:0.001641,0121:0.002219)0.481405:0.001139)0.922387:0.000699)0.995563:0.000238,(0176:0.001720,0022:0.001731)1.000000:0.001179)0.530476:0.000509,0166:0.000812)0.714398:0.000988)0.663217:0.000813,((0015:0.003271,0084:0.003450)0.993189:0.000385,(0125:0.002951,0030:0.002523)0.999517:0.001311)0.999046:0.001242)0.987243:0.000292,(((0156:0.002264,0152:0.000108)0.799351:0.004992,(0066:0.005030,0032:0.003092)0.501974:0.000959)0.840910:0.001364,0042:0.002538)0.999907:0.001394)1.000000:0.000791)1.000000:0.001333,(((0100:0.000692,0025:0.000849)0.988171:0.000188,0104:0.001023)1.000000:0.001904,0079:0.003819)0.424313:0.000791)0.638048:0.001298,0154:0.041527)0.398817:0.000190,((((0175:0.000306,0118:0.000626)1.000000:0.002633,((0148:0.004162,((((((((0086:0.001672,0035:0.001701)0.866223:0.000103,0085:0.000793)0.924941:0.000750,0043:0.001678)0.999439:0.000316,(0054:0.001729,0108:0.001371)0.999658:0.000512)1.000000:0.000355,(0049:0.001504,0017:0.001549)1.000000:0.001159)0.745752:0.001178,0161:0.003015)0.833033:0.001098,0162:0.005303)0.859073:0.001106,0182:0.003297)0.519563:0.000201)0.451509:0.000117,0065:0.002941)0.243742:0.001064)0.838202:0.000141,(0107:0.003409,0039:0.003121)0.968679:0.000820)0.601625:0.003864,0151:0.008490)0.408554:0.000299)0.442116:0.000738,2053:0.007270)0.664370:0.000550,(((0101:0.002529,(0173:0.001014,0018:0.000793)0.971393:0.001023)0.999999:0.001006,0058:0.003910)0.996799:0.000388,(0062:0.002280,(0023:0.002102,0047:0.001917)0.846628:0.000284)1.000000:0.001253)0.927718:0.000786)0.847354:0.000449,0144:0.004764)0.518902:0.000277,0117:0.005346)0.761111:0.000418,(0089:0.001617,0024:0.001890)1.000000:0.001129)0.697174:0.000354)0.657375:0.000095):0.001451,0001:0.001451);
"""  # Replace with full Newick string

# Step 1: Read the map file and create a dictionary
mapping = {}
with open(map_file_path, 'r') as map_file:
    for line in map_file:
        parts = line.strip().split(maxsplit=2)  # Split into number, name, and other data
        if len(parts) >= 2:  # Ensure valid line format
            number, name = parts[0], parts[1]
            mapping[number] = name

# Step 2: Replace numbers in the Newick string
def replace_numbers_in_string(newick_string, mapping):
    return re.sub(r'\b(\d+)\b', lambda match: mapping.get(match.group(1), match.group(0)), newick_string)

# Perform replacements
updated_newick = replace_numbers_in_string(newick_string, mapping)

# Output the updated Newick string
print("Updated Newick String:")
print(updated_newick)
