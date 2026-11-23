# run_all.mpl -- 論文の表をすべて作り直す
#   cd experiments ; maple -q run_all.mpl
# 出力は tables/*.tex（論文から \input される）
read "exp_examples.mpl":
read "exp_illcond.mpl":
read "exp_dense.mpl":
read "exp_random.mpl":
printf("all tables written to tables/\n");
quit;
