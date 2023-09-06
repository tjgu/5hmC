onevar = "amylsqrt"
print(onevar)
b = 10
b = as.numeric(b)
print(b)

j = 9

library(glmnet)

dir = "/orange/zhao/share/5hmC_final/"
dir2 = "./results/residuals/"
dir3 = "./results/residuals/"
dir_o = "results/residuals/"

cl00 = read.table(paste(dir,"Hydro_dataset_655_basic_10-11-2020.Neat1060.txt", sep=""), sep="\t", header=T)
hmc0 = read.table(paste(dir2, "residual_TMM_", onevar, "_scaled_HFv", j, "_others.txt", sep=""), header=T)
chrX = grep("C23", rownames(hmc0))
chrY = grep("C24", rownames(hmc0))
hmc = hmc0[-c(chrX, chrY), ]
colnames(hmc) = sub("X", "", colnames(hmc))

cl = cl00[match(colnames(hmc), cl00$ColName_Hydro), ]

ros_ids = grep("ROS", cl$study)
map_ids = grep("MAP", cl$study)
cl_ros = cl[grep("ROS", cl$study), ]
cl_map = cl[grep("MAP", cl$study), ]

x = apply(hmc, 2, as.numeric)
rownames(x) = rownames(hmc)
x_ros = x[, ros_ids]
x_map = x[, map_ids]

## for map samples
for(k in ((b-1)*10+1):(b*10)){
	a = 0.01
	covs = t(x_map)
	if(onevar == "ad_reagan"){
		one_cl = as.factor(cl[map_ids, which(colnames(cl) == onevar)])
		cvfit <- cv.glmnet(covs, one_cl, alpha=a, nfolds=10, type.measure="mse", family="binomial")
	}else{
		one_cl = as.numeric(cl[map_ids, which(colnames(cl) == onevar)])
		cvfit <- cv.glmnet(covs, one_cl, alpha=a, nfolds=10, type.measure="mse", family="gaussian")
	}
	coef.one = as.matrix(coef(cvfit, s="lambda.min"))
	coef.5hmc = coef.one[2:nrow(coef.one), 1]
	coef.5hmc.sel = coef.5hmc[which(coef.5hmc != 0)]
	
	write.table(coef.5hmc.sel, paste(dir_o, onevar, "_map_resid_v", j,"_a", a, "_k",k, "_5hmc_glmnet_others_rmXY_sel.txt", sep=""), sep="\t", row.names=T, quote=F)
}
