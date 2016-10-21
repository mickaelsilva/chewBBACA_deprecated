#!/usr/bin/python
import CheckCDS
import alleleSizeStats
import os
import argparse
import json
from operator import itemgetter
import HTSeq


def main():
	
	parser = argparse.ArgumentParser(description="This program analyses a set of gene files, analyzing the alleles CDS and the length of the alleles per gene")
	parser.add_argument('-i', nargs='?', type=str, help='list genes, directory or .txt file with the full path', required=True)
	parser.add_argument('-p', nargs='?', type=bool, help='One bad allele still makes gene conserved', required=False)
	parser.add_argument('--log', dest='logScale', action='store_true')
	parser.add_argument('-l', nargs='?', type=str, help='name/location main html file', required=True)
	parser.add_argument('-ta', nargs='?', type=int, help='ncbi translation table', required=True)
	parser.add_argument('-t', nargs='?', type=float, help='Threshold', required=False)
	parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu to use', required=True)
	parser.add_argument('-s', nargs='?', type=int, help='Threshold', required=False)
	parser.set_defaults(logScale=False)
	
	args=parser.parse_args()
	genes = args.i
	transTable = args.ta
	logScale=args.logScale
	htmlFile=args.l
	outputpath=os.path.dirname(htmlFile)
	cpuToUse=args.cpu
	
	try:
		threshold=float(args.t)
	except:
		threshold=0.05
		pass
	try:
		OneBadGeneNotConserved=bool(args.p)
	except:
		OneBadGeneNotConserved=False
		pass
	try:
		splited=int(args.s)
	except:
		splited=False
		pass
	
	try:
		f=open( genes, 'r')
		f.close()
	except IOError:
		listbasename=os.path.basename(os.path.normpath(genes))
		
		with open("listGenes"+listbasename+".txt", "wb") as f:
			for gene in os.listdir(genes):
				try:
					genepath=os.path.join(genes,gene)
					gene_fp2 = HTSeq.FastaReader(genepath)
					for allele in gene_fp2:
						break
					f.write( genepath+"\n")
				except Exception as e:
					print e
					pass
					
				
		genes="listGenes"+listbasename+".txt"
		
	
	genebasename=str(os.path.basename(genes))
	genebasename=genebasename.split(".")
	genebasename.pop()
	genebasename=".".join(genebasename)
	#genebasename=genebasename[0]
	
		
	notConservedgenes,totalgenes,genesWOneAllele,boxplot,histplot,allelenumberplot=alleleSizeStats.getStats(genes,threshold,OneBadGeneNotConserved,True,logScale,outputpath,splited)
	
	#boxplot=str(json.dumps(boxplot))
	histplot=str(json.dumps(histplot))
	allelenumberplot=str(json.dumps(allelenumberplot))


	statsPerGene=CheckCDS.analyzeCDS(genes,transTable,True,outputpath,cpuToUse)
	
	# stats values are ordered in a list allelesNotMultiple3,listStopcodonsInside,listnotStartCodon,numberOfAlleles

	
	with open(htmlFile, "wb") as f:
		f.write("<!DOCTYPE html>\n<html>\n<head><script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js'></script><script type='text/javascript' src='https://mpld3.github.io/js/mpld3.v0.2.js'></script>\n")
		f.write("<script type='text/javascript' src='https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js'></script>")
		f.write("""<!-- Latest compiled and minified JavaScript -->
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>""")
		f.write("""<style type="text/css">
      body {
        padding-top: 10px;
        padding-bottom: 60px;
        padding-left: 20px;
      }

      /* Custom container */
      .container {
        margin: 0 auto;
        max-width: 1000px;
      }
      .container > hr {
        margin: 60px 0;
      }

      /* Main marketing message and sign up button */
      .jumbotron {
        margin: 80px 0;
        text-align: center;
      }
      .jumbotron h1 {
        font-size: 100px;
        line-height: 1;
      }
      .jumbotron .lead {
        font-size: 24px;
        line-height: 1.25;
      }
      .jumbotron .btn {
        font-size: 21px;
        padding: 14px 24px;
      }

    </style>
		<!-- Latest compiled and minified CSS -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">
<script>
function openList1(documentid) {
    var list = document.getElementById(documentid);

    if (list.style.display == "none"){
        list.style.display = "block";
    }else{
        list.style.display = "none";
    }
}
</script>""")
		
		f.write("""<script type='text/javascript'>
    mpld3.register_plugin("clickinfo", ClickInfo);
    ClickInfo.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo.prototype.constructor = ClickInfo;
    ClickInfo.prototype.requiredProps = ["id"];
    ClickInfo.prototype.defaultProps = {labels:null}
    function ClickInfo(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    
    ClickInfo.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        var labels = this.props.labels;
        obj.elements().on("mousedown",function(d, i){ 
                            window.open(labels, '_blank')});
    }
    </script>""")
    
		f.write("""<script type='text/javascript'>
    mpld3.register_plugin("clickinfo2", ClickInfo2);
    ClickInfo2.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo2.prototype.constructor = ClickInfo2;
    ClickInfo2.prototype.requiredProps = ["id"];
    ClickInfo2.prototype.defaultProps = {labels:null}
    function ClickInfo2(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    ClickInfo2.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        labels = this.props.labels;
        obj.elements().on("mousedown",
                          function(d, i){ 
                            window.open(labels[i], '_blank')});
    }
    </script>""")
		
		f.write("""<style type="text/css">
		ul {
    /*min-height: 300px;*/
    -webkit-column-count: 4;
       -moz-column-count: 4;
            column-count: 4; /*4 is just placeholder -- can be anything*/
}
li {
    display: table;
    padding-bottom: 20px; 
    margin-right: 30px;
}
li a {
    color: rgb(0, 162, 232);
}

.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-vn4c{background-color:#D2E4FC;text-align:center}
.tg .tg-qpvr{background-color:#009901;vertical-align:top}
.tg .tg-14d4{background-color:#91df0a;text-align:center;vertical-align:top}

</style>\n</head>\n<body>""")
		
		f.write("""<div class="jumbotron">
  <h1>My Analyzed wg/cg MLST Schema - Rate My Schema</h1>
  <p>Explore the analysis by clicking the analysis buttons</p></div>""")
		

		f.write("<h2>Allele size analysis using a mode +/- "+str(threshold)+"</h2>")
		f.write("<p> Genes are considered not conserved if >1 allele are outside the mode +/-0.05 size. Genes with only 1 allele outside the threshold are considered conserved</p>\n")
		f.write("<h3>"+str(totalgenes)+" total loci</h3>")
		f.write("\n<h3>"+str(len(notConservedgenes))+" loci with high length variability</h3>")
		f.write("""<button onclick = "openList1('ollist1')">Show/hide list</button><ol id='ollist1' style='display: none;'>""")
		
		for elem in notConservedgenes:
			f.write("<li><a href = '"+str(elem)+"' target='_blank'>"+(os.path.basename(elem)).replace(".html","")+"</a></li>")
			
		f.write("</ol>\n<h3>"+str(len(genesWOneAllele))+" loci with only one allele</h3>\n")
		
		f.write("""<button onclick = "openList1('ollist2')">Show/hide list</button><ol id='ollist2' style='display: none;'>""")
		
		for elem in genesWOneAllele:
			f.write("<li><a href = '"+str(elem)+"' target='_blank'>"+(os.path.basename(elem)).replace(".html","")+"</a></li>\n")
		
		f.write("</ol>\n")
		
		f.write("""<div class='container'>
					<div class='row'>
						<div class="col-sm-3">
						<button id='button2' class="btn btn-success btn-block active">Allele length analysis</button>
						</div>
						<div class="col-sm-3">
						<button id='button3' class="btn btn-success btn-block active">Allele number analysis</button>
						</div>
						<div class="col-sm-3">
						<button id='button1' class="btn btn-success btn-block active">Loci size variation analysis</button>
						</div>
						<div class="col-sm-3">
						<button id='button4' class="btn btn-success btn-block active">CDS analysis</button>
						</div>
					</div>
				</div>""")
		
		f.write("""<div id="fig03" style="display:none"><h2>Distribution of number of alleles per gene by mode/mean/median</h2><div id="plot1"></div></div><div id="fig01" style="display:none">
		<h2>Size boxplot for all loci</h2><p>Box plot for each locus on a descending order of the median allele sizes</p>
		<p>Use the zoom button and hover the mouse over a box/median to see the gene name</p>
		<p>-->Blue line represent the median, maximum and minimum</p><p>-->Red line represent the mean</p>
		<button id='buttonbackward' > < </button>
		<button id='buttonforward' > > </button>
		""")
		
		i=0
		for elem in boxplot:
			f.write("""<div id="figbox"""+str(i)+"""" style="display:none"></div>""")
			i+=1
		f.write("""</div><div id="fig02" style="display:none"><h2>Distribution of allele mode sizes per gene</h2></div>""")
		
		f.write("""<script type="text/javascript">
					var boxplotPage=0;
					$("#buttonforward").click(function(){
					  boxplotPage=boxplotPage+1;
					  if ($("#figbox"+boxplotPage).length){
						  if ($("#figbox"+boxplotPage).children().length < 1){
							$.getScript("jsonbox"+boxplotPage+".js");}
						  $("#figbox"+boxplotPage).css({"display":"block"});
						  $("#figbox"+(boxplotPage-1)).css({"display":"none"});
					   }
					   else{
					    boxplotPage=boxplotPage-1;
					    }
					$('html, body').animate({scrollTop:$(document).height()}, 1);
					return false;
					}); 
					</script>""")
		
		f.write("""<script type="text/javascript">
					$("#buttonbackward").click(function(){
					  boxplotPage=boxplotPage-1;
					  if ($("#figbox"+boxplotPage).length){
						  if ($("#figbox"+boxplotPage).children().length < 1){
							$.getScript("jsonbox"+boxplotPage+".js");}
						  $("#figbox"+boxplotPage).css({"display":"block"});
						  $("#figbox"+(boxplotPage+1)).css({"display":"none"});
					   }
					   else{
					    boxplotPage=boxplotPage+1;
					    }
					$('html, body').animate({scrollTop:$(document).height()}, 1);
					return false;
					}); 
					</script>""")
		
		f.write("""<script type="text/javascript">
					$("#button2").click(function(){
					  if ($("#fig02 .mpld3-figure").length < 1){
						$.getScript("json2.js");}
					  $("#fig02").css({"display":"block"});
					  $("#fig01").css({"display":"none"});
					  $("#fig03").css({"display":"none"});
					  $("#fig04").css({"display":"none"});
					}); 
					</script>""")
		
		f.write("""<script type="text/javascript">
					$("#button3").click(function(){
					if ($("#fig03 .mpld3-figure").length < 1){
						$.getScript("json3.js");}
					  $("#fig03").css({"display":"block"});
					  $("#fig02").css({"display":"none"});
					  $("#fig01").css({"display":"none"});
					  $("#fig04").css({"display":"none"});
					}); 
					</script>""")
		
		f.write("""<script type="text/javascript">
					$("#button1").click(function(){
					if ($("#fig01 .mpld3-figure").length < 1){
						$.getScript("jsonbox0.js");}
					  $("#fig01").css({"display":"block"});
					  $("#figbox0").css({"display":"block"});
					  $("#fig02").css({"display":"none"});
					  $("#fig03").css({"display":"none"});
					  $("#fig04").css({"display":"none"});
					}); 
					</script>""")
		f.write("""<script type="text/javascript">
					$("#button4").click(function(){
					  $("#fig04").css({"display":"block"});
					  $("#fig02").css({"display":"none"});
					  $("#fig03").css({"display":"none"});
					  $("#fig01").css({"display":"none"});
					}); 
					</script>""")
		
		
		f.write("""<div id="fig04" style="display:none"><title>Schema Validation Results</title>\n<h1>Allele CDS analysis results</h1>\n<p>Summary table of the alleles with issues per gene using the <a href='http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11'>NCBI translation table 11</a></p><p>Click on the gene name to open the fasta file</p><p>click on the boxes with the % to get the index of the alleles with issues</p><table class="tg">
  <tr>
    <th class="tg-031e">Gene</th>
    <th class="tg-031e">Number alleles not multiple of 3</th>
    <th class="tg-031e">Number alleles w/ >1 stop codons</th>
    <th class="tg-031e">Number alleles wo/ Start/Stop Codon</th>
    <th class="tg-qpvr" .background-color='#59b300'>Number of alleles (% alleles w/ issues) </th>
  </tr>""")
		ordered=[]
		for key, value in statsPerGene.iteritems():
			aux=[]
			numberMultip=float(len(value[0]))
			numberStop=float(len(value[1]))
			numberStart=float(len(value[2]))
			total=float(value[3])
			totalpercent=((numberMultip+numberStart+numberStop)/total)*100
			aux.append(key)
			aux.append(totalpercent)
			ordered.append(aux)
		ordered=sorted(ordered, key=itemgetter(-1))
		ordered.reverse()
		
		newlist=[]	
		for item in ordered:	
			
			aux=[]
			aux.append(item[0])
			
			
			
			value=statsPerGene[item[0]]
			numberMultip=float(len(value[0]))
			numberStop=float(len(value[1]))
			numberStart=float(len(value[2]))
			total=float(value[3])
			
			aux.append(value)
			newlist.append(aux)
			name=os.path.basename(str(item[0]))
			name=name.split(".")
			name=name[0]
			if (numberMultip>0 or numberStop>0 or numberStart>0):
				f.write("<tr id="+str(item[0])+""">\n<td class='tg-vn4c' onclick="window.open('"""+str(item[0])+"""')" style='cursor:pointer'>"""+name+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(numberMultip))+" ("+str('{0:.2f}'.format((numberMultip/total)*100))+"%)"+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(numberStop))+" ("+str('{0:.2f}'.format((numberStop/total)*100))+"%)"+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(numberStart))+" ("+str('{0:.2f}'.format((numberStart/total)*100))+"%)"+"</td>\n<td class='tg-14d4' onclick='a(this);'>"+str(int(total))+" ("+str('{0:.2f}'.format(((numberMultip+numberStart+numberStop)/total)*100))+"%)"+"</td>\n</tr>")
		
		
		f.write("</table>")
		f.write("""<div id='AllelesWissues'></div><button onclick="$('#AllelesWissues').empty();">clean</button></div>""")
		f.write("""\n<script type='text/javascript'>function a(element) {
	var id = $(element).closest("tr").attr("id");
	var badalleles=[];
	for (i = 0; i < alleles.length; i++) { 
		if((alleles[i])[0]==id){
			badalleles=(alleles[i])[1];
			break;
			}
		}
	var notmulti=(badalleles[0]).join('; ');
	var stopcodon=(badalleles[1]).join('; ');
	var startcodon=(badalleles[2]).join('; ');
	var name=(id.split("/")).slice(-1)[0]
	name=(name.split("."))[0]
	$('#AllelesWissues').append('<h2> Gene: '+name+'</h2>');
	$('#AllelesWissues').append('<p> Alleles not multiple of 3: '+notmulti+'</p><p> Alleles with >1 stop codon: '+stopcodon+'</p><p>Alleles without start codon: '+startcodon+'</p>');

	$('html,body').animate({
        scrollTop: $('#AllelesWissues').offset().top},'slow');
	
	}
	</script>""")
		
		f.write("\n<script type='text/javascript'>var alleles="+json.dumps(newlist)+"</script>")
		f.write("</body>\n</html>")
	
	i=0
	for elem in boxplot:
		
		boxplotElem=str(json.dumps(elem))
		filename="jsonbox"+str(i)+(".js")
		with open((os.path.join(outputpath,filename)), "wb") as f:
			
			f.write("var jsonbox"+str(i)+" ="+str(boxplotElem)+";mpld3.draw_figure('figbox"+str(i)+"', jsonbox"+str(i)+");")
		i+=1
		
	with open((os.path.join(outputpath,"json2.js")), "wb") as f:
		
		f.write("var json02 ="+str(histplot)+";mpld3.draw_figure('fig02', json02);")
	
	with open((os.path.join(outputpath,"json3.js")), "wb") as f:
		
		f.write("var json03 ="+str(allelenumberplot)+";mpld3.draw_figure('plot1', json03);")


	try:
		os.remove("listGenes"+listbasename+".txt")
	except:
		pass
	
if __name__ == "__main__":
    main()

