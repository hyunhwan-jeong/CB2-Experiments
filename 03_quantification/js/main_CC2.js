var lineByLine = require('n-readlines');

var liner_for_library = new lineByLine(process.argv[2]);

var library_map = new Map();
var count_map = new Map();
var lib_seq_len = 0;
while(line = liner_for_library.next()) {
    var line = line.toString().split('\t');
    var num = 0;

    for(var i=0;i<line[1].length;++i) {
        num *= 4;
        num += (line[1].charCodeAt(i)>>1)&3;
    }
    lib_seq_len = line[1].length;
    //console.log(line[1] + " " + num);
    library_map.set(line[0], num);
    count_map.set(num, 0);
}

var liner_for_fastq = new lineByLine(process.argv[3]);

var line_number = 0;
var num_proceed_lines = 0;
var num_hits = 0;
var mod = 4**lib_seq_len;

while(line = liner_for_fastq.next()) {
    if(line_number++%4==1) {
        line = line.toString();
        if(++num_proceed_lines%1000000==0) {
            console.info("Processing "+ num_proceed_lines+"th lines...");
        }

        var j = 0; 
        var num = 0;
        var len = line.length;
        for(var i = 0 ; i < len ; ++i) {
            if(line[i]=='N') {
                j = 0; num = 0;
                continue;
            }

            var p = (line.charCodeAt(i)>>1)&3;
            num *= 4;
            num += p;
            num = num % mod;
            if( ++j == lib_seq_len ) {
                if( count_map.has(num) ) {
                    ++num_hits;
                    count_map.set(num, count_map.get(num)+1);
                    break;
                }
                j--;
            }

        }
    }
}

console.log(num_hits+"/"+num_proceed_lines);
for(var gRNA of library_map.keys()) {
    var seq = library_map.get(gRNA);
    var count = count_map.get(seq);
    console.log(gRNA+"\t"+count);
}
