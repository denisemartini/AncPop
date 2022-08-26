#produce VCFs for the make_fasta script
#use the command format:
#	filter_cover_for_vcfs.R [vcf name] <min coverage>
#eg: filter_cover_for_vcfs.R harrier_coverage_gatk_table Harrier.test-nohead.vcf 3

#specify arguments from the command line
args<-commandArgs(TRUE)
if(length(args)<2){
	print("include a coverage table, vcf without headers and coverage minimum")
	stopifnot(length(args)>1)
} else if(length(args)==2){
	print("Producing coverage with default minimum coverage of 1")
	coverage.table<-args[1]
	vcf.nohead<-args[2]
	min.cov<-1
} else if(length(args)==3){
	coverage.table<-args[1]
	vcf.nohead<-args[2]
	min.cov<-as.numeric(args[3])
} else {
	print("Too many arguments")
	stopifnot(length(args)<4)
}

coverage<-read.table(coverage.table,header=F,stringsAsFactors=F)

gens<-read.table(vcf.nohead,header=F)

gen.coors<-paste0(gens[,1],":",gens[,2])
cov.coors<-paste0(coverage[,1],":",coverage[,2])

#cov.coors<-coverage[,1]

#gen.coors[!(gen.coors%in%cov.coors)]

#cov.coors[!(cov.coors%in%gen.coors)]

cov.tab.min<-coverage[coverage[,3]<min.cov,]

cov.tab.min.coor<-paste0(cov.tab.min[,1],":",cov.tab.min[,2])

gens.min.cov<-gens[gen.coors%in%cov.tab.min.coor,]

write.table(gens.min.cov,paste0("gens.",min.cov,".vcf"),quote=F,col.names=F,row.names=F,sep="\t")
