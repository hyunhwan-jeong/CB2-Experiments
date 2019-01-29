var lineByLine = require('n-readlines');

var liner_for_library = new lineByLine(process.argv[2]);

var library_map = new Map();
var count_map = new Map();
var lib_seq_len = 0;
while(line = liner_for_library.next()) {
    var line = line.toString().split('\t');
    lib_seq_len = line[1].length;
    library_map.set(line[0], line[1]);
    count_map.set(line[1], 0);
}

var liner_for_fastq = new lineByLine(process.argv[3]);

var line_number = 0;
var num_proceed_lines = 0;
var num_hits = 0;
var max_sample_line = 100000;

while(line = liner_for_fastq.next()) {
    if(line_number++%4==1) {
        line = line.toString();
        if(++num_proceed_lines%1000000==0) {
            console.info("Processing "+ num_proceed_lines+"th lines...");
        }
        for(var i = 0 ; i+lib_seq_len <= line.length ; ++i ) {
            var sub = line.substring(i,i+lib_seq_len);
            if(count_map.has(sub)) {
                count_map.set(sub, count_map.get(sub)+1);
                ++num_hits;
                break;
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
