 #!/bin/sh

##################################################################
#   Script for running a MPI job @ GCOE saver
#
# ___ Notice ___
#   (1) Execusion command -> sh go.sh &
#   (2) Attribute of this file must be '744'
#   (3) 一回実行すると、もう一度ログインし直さないと変数はリセットされない
#   (4) 実行ファイルは広域計算用の設定でMakefileすること
##################################################################

#Deleting previous output files
   echo "Deleting previous output files"
   cd result
   rm      *_out.txt          #Delete previous result files 1
   rm     *_out2.txt          #Delete previous result files 2
#   rm     *spnin.txt         #(リスタート時にのみ意味あり) 古いスピンインファイルを削除
#   rename spnout spnin *.txt #(リスタート時にのみ意味あり) スピンアウトファイルをスピンインファイルに名前変更
   cd ../
   
#Submit MPI program, and wait for complete
   qsub go.bat
   echo "Waiting for MPI comuptation"
   date
   while [[ `qstat | grep hsato | wc -l` -ge 1  ]]; do
     sleep 5
   done
   date
   
#結果ファイルの変換プログラムを投入し、これが終わるまで待つ
#   cd code_visualize
#   qsub go.bat
#   cd ../
#   echo "Converting result files"
#   
#   while [[ `qstat | grep hsato | wc -l` -ge 1  ]]; do
#     sleep 3
#   done
   
#Visualizing result files with R script
   echo "Visualizing result files"
   R --vanilla --slave --quiet < visualization.R
   
#Delete job record files
   rm *.o[0-9]*
   
#message
   echo "Normal End"
