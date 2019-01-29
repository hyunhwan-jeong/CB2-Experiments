var lineByLine = require('n-readlines');

var liner_for_library = new lineByLine("library/nbt3536-S3.tsv");

var library_map = new Map();
var count_map = new Map();
while(line = liner_for_library.next()) {
    var line = line.toString().split('\t');
    var num = 0;

    for(var i=0;i<line[1].length;++i) {
        num *= 4;
        num += (line[1].charCodeAt(i)>>1)&3;
    }
    //console.log(line[1] + " " + num);
    library_map.set(line[0], num);
    count_map.set(num, 0);
}

var liner_for_fastq = new lineByLine("SRR3341899.fastq");

var line_number = 0;
var num_proceed_lines = 0;
var num_hits = 0;
while(line = liner_for_fastq.next()) {
    if(line_number++%4==1) {
        line = line.toString();
        if(++num_proceed_lines%100000==0) {
            console.info("Processing "+ num_proceed_lines+"th lines...");
        }

        var j = 0; 
        var num = 0;
        for(var i = 0 ; i < line.length ; ++i) {
            if(line[i]=='N') {
                j = 0; num = 0;
                continue;
            }

            var p = (line.charCodeAt(i)>>1)&3;
            num *= 4;
            num += p;
            num = num % 1099511627776;
            if( ++j == 20 ) {
                if( count_map.has(num) ) {
                    count_map.set(num, count_map.get(num)+1);
                    break;
                }
                j--;
            }

        }
    }
}

for(var gRNA of library_map.keys()) {
    var seq = library_map.get(gRNA);
    var count = count_map.get(seq);
    console.log(gRNA+"\t"+count);
}
