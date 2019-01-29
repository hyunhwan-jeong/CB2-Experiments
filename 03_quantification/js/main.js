var lineByLine = require('n-readlines');

var liner_for_library = new lineByLine("library/nbt3536-S3.tsv");

var library_map = new Map();
var count_map = new Map();
while(line = liner_for_library.next()) {
    var line = line.toString().split('\t');
    library_map.set(line[0], line[1]);
    count_map.set(line[1], 0);
}

var liner_for_fastq = new lineByLine("SRR3341899.fastq");

var line_number = 0;
var num_proceed_lines = 0;
while(line = liner_for_fastq.next()) {
    if(line_number++%4==1) {
        line = line.toString();
        if(++num_proceed_lines%100000==0) {
            console.info("Processing "+ num_proceed_lines+"th lines...");
        }

        for(var i = 0 ; i+20 <= line.length ; ++i) {
            var sub_str = line.substring(i,i+20);
            if(count_map.has(sub_str)) {
                count_map.set(sub_str, count_map.get(sub_str)+1);
                break;
            }
        }
    }
}

for(var gRNA of library_map.keys()) {
    var seq = library_map.get(gRNA);
    var count = count_map.get(seq);
    console.log(gRNA+"\t"+count);
}


