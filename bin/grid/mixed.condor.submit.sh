######################################################################
# 				MixedEvent
######################################################################

InitialDir = /home/jdb12/work/FemtoDstPairAna/bin/
Executable = /home/jdb12/work/FemtoDstPairAna/bin/pairAna.app
Arguments  = /home/jdb12/work/FemtoDstPairAna/bin/config/MixedEvent_invPID.xml --jobIndex=$(Process)

GetEnv     = True

Queue 25