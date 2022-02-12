import os, commands, re
import sqlite3
# outFolder = '/storage/af/group/rdst_analysis/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples_B2DstMu_wOC'
#
# jobsId = '875614'
#
#
# cmd = 'condor_q {} -hold'.format(jobsId)
# status, output = commands.getstatusoutput(cmd)
#
# failedIds = []
# for l in output.split('\n'):
#     if re.match('[0-9]+\.[0-9]+[ ]+ocerri', l):
#         failedIds.append(l.split(' ')[0].split('.')[1])
#
# print 'Total failed jobs', len(failedIds)
#
# additionalFailingIds = []
# # 712, 747, 786, 763, 847, 966, 971,]
#
# for id in failedIds + [str(x) for x  in additionalFailingIds]:
#     # if int(id) == 1000:
#     #     break
#     print id
#     cmd = 'rm {}/out_CAND_{}.root'.format(outFolder, id)
#     # print cmd
#     os.system(cmd)
#     cmd = 'rm {}/out/job_{}_{}.*'.format(outFolder, id, jobsId)
#     # print cmd
#     os.system(cmd)

home = os.path.expanduser("~")
db = os.path.join(home,'state.db')

conn = sqlite3.connect(db)

c = conn.cursor()
results = c.execute('SELECT id, log_file, output_file FROM ntuplizer_jobs WHERE state="FAILED"')

nRemoved_root = 0
nRemoved_log = 0
for row in results.fetchall():
    id, log_file, output_file = row
    if os.path.isfile(output_file):
        print 'Removing output of failed job', id
        os.remove(output_file)
        nRemoved_root += 1
    if os.path.isfile(log_file):
        print 'Removing log of failed job', id
        os.remove(log_file.replace('.log', '.err'))
        os.remove(log_file.replace('.log', '.out'))
        os.remove(log_file)
        nRemoved_log += 1

print 'Total removed root files:', nRemoved_root
print 'Total removed log files:', nRemoved_log
