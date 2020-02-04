# You have to install the library function first
install.packages("rvest")
install.packages("enrichR")

# Import the library
library(rvest) #For web crawler

# function of removing the redundant ones
Delete_Double=function(input){
	n=length(input)/2
	output=rep(NA,n)
	for (i in 1:n){
		output[i]=input[1+(i-1)*2]
	}
	return(output)
}

#microRNA=readline(prompt="Enter Name of microRNA:")

microRNA="hsa-miR-144-3p"; #Name of microRNA

url=paste("http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&miRNAs%5B%5D=","&genes%5B%5D=&species%5B%5D=1&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1",sep=microRNA) #Link of Tarbase website



doc=read_html(paste("http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&miRNAs%5B%5D=","&genes%5B%5D=&species%5B%5D=1&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1",sep=microRNA)) #Crawl the website

# Find the names of target genes
temp=doc%>%html_nodes(".first-level-block-bold")%>%html_text()
temp2=gsub(" \n","",temp) #substrings
results=Delete_Double(temp2) #Removing the redundant ones

write.table(results, file = "/Users/user/Desktop/TarBase_Genes.csv",sep = ",", quote = FALSE, na = "NA") #Outputting the Gene list to the Desktop

#################### Using Enrichr(GO) in R ####################
library(enrichR) #For GO analysis
lists=enrichr(results, databases = "GO_Biological_Process_2018") #Using GO_Biological_Process_2018 database as example
#Link of databases->http://amp.pharm.mssm.edu/Enrichr/#stats
df <- ldply (lists, data.frame)
df = =df[1:20,] # Get top 20 of GO list
