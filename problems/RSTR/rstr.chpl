module RSTR {
  use Math;
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  proc string_log_prob(sequence, gcprob) {
    var nbases = ["A" => 0, "C" => 0, "G" => 0, "T" => 0];
    for base in sequence do
        nbases[base] += 1;

    const loggprob = log(gcprob/2.0);
    const logaprob = log((1-gcprob)/2.0);

    var logprob = nbases["G"]*loggprob + nbases["C"]*loggprob;
    logprob += nbases["A"]*logaprob + nbases["T"]*logaprob;

    return logprob;
  }

  proc main() {
    try! {
      var channel = input_channel(infile, defaultfilename);
      var seq, values : string;

      channel.readline(values);
      const items = values.split();
      var n : int = items[1] : int;
      var gcprob : real = items[2] : real;
      
      channel.readline(seq);
      seq.strip();

      const sequence_non_prob = -expm1(string_log_prob(seq, gcprob));
      writeln(1.0 - sequence_non_prob**n);
    }
  }
}
