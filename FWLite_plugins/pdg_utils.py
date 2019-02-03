pdgId2Name_dic = {}

pdgId2Name_dic[11] = 'e-'
pdgId2Name_dic[-11] = 'e+'
pdgId2Name_dic[12] = 'nu_e'
pdgId2Name_dic[-12] = 'anti-nu_e'
pdgId2Name_dic[13] = 'mu-'
pdgId2Name_dic[-13] = 'mu+'
pdgId2Name_dic[14] = 'nu_u'
pdgId2Name_dic[-14] = 'anti-nu_u'
pdgId2Name_dic[15] = 'tau-'
pdgId2Name_dic[-15] = 'tau+'
pdgId2Name_dic[16] = 'nu_t'
pdgId2Name_dic[-16] = 'anti-nu_t'

pdgId2Name_dic[21] = 'g'
pdgId2Name_dic[22] = 'gamma'

pdgId2Name_dic[111] = 'pi0'
pdgId2Name_dic[211] = 'pi+'
pdgId2Name_dic[-211] = 'pi-'

pdgId2Name_dic[321] = 'K+'
pdgId2Name_dic[-321] = 'K-'

pdgId2Name_dic[411] = 'D+'
pdgId2Name_dic[-411] = 'D-'
pdgId2Name_dic[413] = 'D*+'
pdgId2Name_dic[-413] = 'D*-'
pdgId2Name_dic[421] = 'D0'
pdgId2Name_dic[-421] = 'anti-D0'
pdgId2Name_dic[423] = 'D0*'
pdgId2Name_dic[-423] = 'anti-D0*'

pdgId2Name_dic[511] = 'B0'
pdgId2Name_dic[-511] = 'anti-B0'
pdgId2Name_dic[521] = 'B+'
pdgId2Name_dic[-521] = 'B-'
pdgId2Name_dic[523] = 'B*+'
pdgId2Name_dic[-523] = 'B*-'

def getName(pdgId):
    if pdgId in pdgId2Name_dic.keys():
        return str(pdgId2Name_dic[pdgId])
    else:
        return str(pdgId)

pdgId2MeanLife_dic = {} #ns
pdgId2MeanLife_dic[511] = 1.520e-3
pdgId2MeanLife_dic[521] = 1.638e-3

def getMeanLife(pdgId):
    if abs(pdgId) in pdgId2MeanLife_dic.keys():
        return pdgId2MeanLife_dic[abs(pdgId)]
    else:
        print 'Add Mean life to dic for', abs(pdgId)
        exit()

pdgId2mass = {} #ns
pdgId2mass[511] = 5.2796
pdgId2mass[521] = 5.2793

def getMass(pdgId):
    if abs(pdgId) in pdgId2mass.keys():
        return pdgId2mass[abs(pdgId)]
    else:
        print 'Add Mean life to dic for', abs(pdgId)
        exit()
