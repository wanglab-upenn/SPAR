submoduleDir=$1
ls ${submoduleDir}/*.r | \
  awk 'BEGIN{printf "<html>\n<body>\n";}
     {
       plotname=$0;
       gsub(/^.+ggplot2_/,"",plotname);
       gsub(/.r$/,"",plotname);
       img=("<img src=\"" plotname ".png\" width=\"800\">");
       print img
     }END{printf "</body>\n</html>"}' 
#> "${OUTDIR}/${module}_${step}.html"

