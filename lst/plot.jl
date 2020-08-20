using CSV
using Plots
using LaTeXStrings
record2 = CSV.read("tab/2record.csv")
scatter(
xlabel = "序号",
ylabel = L"v/(\mathrm{m/s})",
label = ["标定速度" "测量速度" "绝对误差"],
xtick = record2[1],
[record2[6], record2[7], record2[8]]
)
savefig("fig/2record.pdf")
record3 = CSV.read("tab/3record.csv")
scatter(
xlabel = "序号",
ylabel = L"r/\mathrm{m}",
label = ["标定距离" "测量距离" "绝对距离"],
xtick = record2[1],
[record3[4], record3[5], record3[6]]
)
savefig("fig/3record.pdf")
