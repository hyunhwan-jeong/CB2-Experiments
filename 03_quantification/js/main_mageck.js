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

var trimkey = new Map();


while(line = liner_for_fastq.next()) {
    if(line_number++%4==1) {
        line = line.toString();
        if(++num_proceed_lines%100000==0) {
            console.info("Processing "+ num_proceed_lines+"th lines...");
        }
        if(num_proceed_lines>max_sample_line) break;
        for(var i = 0 ; i+lib_seq_len <= line.length ; ++i) {
            var sub = line.substring(i,i+lib_seq_len);
            if(count_map.has(sub)) {
                if(!trimkey.has(i)) trimkey.set(i,0);
                trimkey.set(i, trimkey.get(i)+1);
                break;
            }
        }
    }
}

var total_mapped_reads = [...trimkey.values()].reduce((a,b)=>a+b,0);
var sorted =[...trimkey.entries()].sort(function(a,b) {
    if(a[1]!=b[1]) return b[1]-a[1];
    else return a[0]-b[0];
});

var cands = [];
var last_frac = 0.01;
var cum_frac = 0;
for(var cand of sorted) {
    var cur_frac = 1.*cand[1]/total_mapped_reads;
    console.log(cand[0]+"\t"+cand[1]+"\t"+cur_frac);
    if(last_frac>cur_frac*3) break;
    cum_frac += cur_frac;
    cands.push(cand[0]);
    if(cum_frac>0.75) break;
    if(cands.length>0 && cur_frac<0.05) break;
    last_frac = cur_frac;
}
console.log("Candidates: "+ cands);


liner_for_fastq = new lineByLine(process.argv[3]);

line_number = 0;
num_proceed_lines = 0;
num_hits = 0;


while(line = liner_for_fastq.next()) {
    if(line_number++%4==1) {
        line = line.toString();
        if(++num_proceed_lines%1000000==0) {
            console.info("Processing "+ num_proceed_lines+"th lines...");
        }
        for(var i of cands) {
            if( i+lib_seq_len < line.length ) {
                var sub = line.substring(i,i+lib_seq_len);
                if(count_map.has(sub)) {
                    count_map.set(sub, count_map.get(sub)+1);
                    ++num_hits;
                    break;
                }
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
