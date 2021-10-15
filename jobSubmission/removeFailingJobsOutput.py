import os, commands, re

outFolder = '/storage/af/group/rdst_analysis/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples_B2DstMu_wOC'

jobsId = '875614'


cmd = 'condor_q {} -hold'.format(jobsId)
status, output = commands.getstatusoutput(cmd)

failedIds = []
for l in output.split('\n'):
    if re.match('[0-9]+\.[0-9]+[ ]+ocerri', l):
        failedIds.append(l.split(' ')[0].split('.')[1])

print 'Total failed jobs', len(failedIds)

additionalFailingIds = []
# 712, 747, 786, 763, 847, 966, 971,]

for id in failedIds + [str(x) for x  in additionalFailingIds]:
    # if int(id) == 1000:
    #     break
    print id
    cmd = 'rm {}/out_CAND_{}.root'.format(outFolder, id)
    # print cmd
    os.system(cmd)
    cmd = 'rm {}/out/job_{}_{}.*'.format(outFolder, id, jobsId)
    # print cmd
    os.system(cmd)
