module LIA {
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const off_per_couple = 2 : uint;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) {
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

    for i in 0..2 do
      for j in 0..2 {
        T0[i, j] = (2-i)*(2-j)/4.0;
        T2[i, j] = i*j/4.0;
      }

    allele_probs[0] = dot( parent2, dot( T0, parent1 ));
    allele_probs[2] = dot( parent2, dot( T2, parent1 ));
    allele_probs[1] = 1.0 - allele_probs[0] - allele_probs[2];
    return allele_probs;
  }

  proc frac_expected_dominant(parent1, parent2) {
    var probs = allele_count_probabilities(parent1, parent2);
    return probs[2] + probs[1];
  }

  proc main() {
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
