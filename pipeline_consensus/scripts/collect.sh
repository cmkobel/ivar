#!/bin/bash


f_pangolin="output/all_pangolin.csv"
f_nextclade="output/all_nextclade.csv"
f_input="output/all_input_lists.tab"
f_coverage="output/all.coverage.tab"


# collect pangolin 
echo "pangolin ..."
echo "sample,taxon,lineage,SH-alrt,UFbootstrap,lineages_version,status,note" \
> $f_pangolin
cat output/*/pangolin/*pangolin.csv | grep -v "taxon,lineage," >> $f_pangolin



# collect nextclade
echo "nextclade ..."
echo "sample;seqName;substitutions;totalMutations;aminoacidChanges;totalAminoacidChanges;insertions;totalInsertions;deletions;totalGaps;missing;totalMissing;nonACGTNs;totalNonACGTNs;alignmentStart;alignmentEnd;alignmentScore;pcrPrimerChanges;totalPcrPrimerChanges;clade;qc.seqName;qc.privateMutations.score;qc.privateMutations.total;qc.privateMutations.excess;qc.privateMutations.cutoff;qc.privateMutations.status;qc.missingData.score;qc.missingData.totalMissing;qc.missingData.missingDataThreshold;qc.missingData.status;qc.snpClusters.score;qc.snpClusters.totalSNPs;qc.snpClusters.clusteredSNPs;qc.snpClusters.status;qc.mixedSites.score;qc.mixedSites.totalMixedSites;qc.mixedSites.mixedSitesThreshold;qc.mixedSites.status;qc.overallScore;qc.overallStatus" \
> $f_nextclade
cat output/*/nextclade/*_nextclade.csv | grep -v "seqName;substitutions;" >> $f_nextclade


# collect input list files
echo "input list ..."
#echo -e "forward\\treverse\\tsample\\tcomment" > $f_input
echo -e "forward\\treverse\\tsample\\tcomment\\tforward_path\\treverse_path\\tbatch\\tfull_name" > $f_input
cat output/*_input_list.tab | grep -vP "forward\\treverse" >> $f_input


# collect coverage list files
echo "coverage ..."
echo -e "sample\\tfilter\\ttarget\\tposition\\tcoverage" > $f_coverage
cat output/*_coverage.tab >> $f_coverage





# Call R and merge the tables together? Or is that too much?
