
cat bedForROHAnalysis_ROH.hom.indiv | awk '{print "XXFID", "XXIID",$3,$4,$5,$6}' | head -n1 > annonymized_bedForROHAnalysis_ROH.hom.indiv

cat bedForROHAnalysis_ROH.hom.indiv | awk '{print "XXFID", "XXIID",$3,$4,$5,$6}' | tail -n+2 | shuf -n200 >> annonymized_bedForROHAnalysis_ROH.hom.indiv

cat bedForROHAnalysis_ROH.hom | sed "s/ * / /g" | awk '{print "XXFID", "XXIID",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | head -n1 > annonymized_bedForROHAnalysis_ROH.hom

cat bedForROHAnalysis_ROH.hom | sed "s/ * / /g" | awk '{print "XXFID", "XXIID",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' |  tail -n+2 | shuf -n2000  >> annonymized_bedForROHAnalysis_ROH.hom

