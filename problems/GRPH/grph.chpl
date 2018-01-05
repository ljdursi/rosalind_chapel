module GRPH {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const overlap = 3;

  record read {
    var seq: string;
    var readlabel: string;
  }

  record array_of_reads {
    var D: domain(1);
    var reads : [D] read;
  }

  // use stdin if filename == defaultfilename
  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
    var sequence = "", seqlabel = "";
    var started = false;
    var line = "";

    while (inchannel.readline(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..).strip();
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var allreads : domain(read);
    
      var kmers : domain(string);
      var prefixes : [kmers] array_of_reads;
  
      for (seqlabel, sequence) in fasta_iterator(inchannel) {
        // populate the nodes
        const curread = new read(sequence, seqlabel);
        allreads.add(curread);
        prefixes[sequence(1..overlap)].reads.push_back(curread);
      }
     
      for read in allreads {
        const l = read.seq.len;
        const suffix = read.seq(l-overlap+1..l);
        if kmers.member(suffix) then
          for neighbor in prefixes[suffix].reads do
            if neighbor.readlabel != read.readlabel then
              writeln(read.readlabel, " ", neighbor.readlabel);
      }
    }
  }
}
