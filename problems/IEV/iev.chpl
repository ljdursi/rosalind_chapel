module IEV {
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const off_per_couple = 2 : uint;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var counts: [1..6] uint;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(counts);
    return counts;
  }

  proc allele_count_probabilities(parent1, parent2) {
    var T0 = Matrix(3, 3);  
    var T2 = Matrix(3, 3);  
    var allele_probs = Vector(3), curprobs = Vector(3);

    for i in 1..3 do
      for j in 1..3 {
        T0[i, j] = (3-i)*(3-j)/4.0;
        T2[i, j] = (i-1)*(j-1)/4.0;
      }

    allele_probs[1] = dot( parent2, dot( T0, parent1 ));
    allele_probs[3] = dot( parent2, dot( T2, parent1 ));
    allele_probs[2] = 1.0 - allele_probs[1] - allele_probs[3];
    return allele_probs;
  }

  proc frac_expected_dominant(parent1, parent2) {
    var probs = allele_count_probabilities(parent1, parent2);
    return probs[3] + probs[2];
  }

  proc main() {
    try! {
      var counts = params_from_file(infile, defaultfilename);
      var AA = Vector([0.0, 0.0, 1.0]);
      var Aa = Vector([0.0, 1.0, 0.0]);
      var aa = Vector([1.0, 0.0, 0.0]);

      var expectation = 0.0;
      expectation += off_per_couple * counts[1] * frac_expected_dominant(AA, AA);
      expectation += off_per_couple * counts[2] * frac_expected_dominant(AA, Aa);
      expectation += off_per_couple * counts[3] * frac_expected_dominant(AA, aa);
      expectation += off_per_couple * counts[4] * frac_expected_dominant(Aa, Aa);
      expectation += off_per_couple * counts[5] * frac_expected_dominant(Aa, aa);
      expectation += off_per_couple * counts[6] * frac_expected_dominant(aa, aa);

      writef("%dr\n", expectation);
    }
  }
}
