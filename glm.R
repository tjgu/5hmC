onevar = "amylsqrt"
j = 9
print(onevar)
print(j)

dir = "/orange/zhao/share/5hmC_final/"
dir2 = "./results/residuals/"
dir_o = "./results/residuals/"

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
res = matrix(NA, nrow(x_map), 4)
rownames(res) = rownames(x_map)
colnames(res) = c("Qval", "Pval", "logOR", "OR")

if(onevar == "ad_reagan"){
        one_cl = as.factor(cl[map_ids, which(colnames(cl) == onevar)])
}else{
        one_cl = as.numeric(cl[map_ids, which(colnames(cl) == onevar)])
}
toRM = which(is.na(one_cl))
if(length(toRM) > 0){
	for(k in 1:nrow(x_map)){
		if(onevar == "ad_reagan"){
			onelg = glm(one_cl[-toRM] ~ x_map[k, -toRM], family=binomial(link='logit'))
		}else{
			onelg = glm(one_cl[-toRM] ~ x_map[k, -toRM])
		}
		onelgsum = summary(onelg)$coefficients
		res[k, 2:3] = c(onelgsum[2,4], onelgsum[2,1])
	}
}else{
	for(k in 1:nrow(x_map)){
		if(onevar == "ad_reagan"){
			onelg = glm(one_cl ~ x_map[k, ], family=binomial(link='logit'))
		}else{
			onelg = glm(one_cl ~ x_map[k, ])
		}
		onelgsum = summary(onelg)$coefficients
		res[k, 2:3] = c(onelgsum[2,4], onelgsum[2,1])
        }
}
res[,4] = exp(res[,3])
qval = p.adjust(res[,2], method="bonferroni")
res[,1] = qval

print(paste("qval <0.05: ", length(which(qval < 0.05))))
print(paste("summary of OR: "))
print(summary(res[,3]))
print(summary(res[which(res[,1] < 0.05), 3]))
write.table(res, paste(dir_o, onevar, "_map_resid_v", j, "_5hmc_others_rmXY.txt", sep=""), sep="\t", row.names=T, quote=F)
print(paste("finish", j))
gc()

