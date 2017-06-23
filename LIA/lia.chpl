module LIA {
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const nalleles = 2;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) {
    var k: uint;
    var n: uint;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(k, n);
    return (k, n);
  }

  proc allele_count_probabilities(ngen, original, other) {
    var T0 = Matrix(3, 3);  
    var T2 = Matrix(3, 3);  
    var allele_probs = Vector(3), curprobs = Vector(3);

    for i in 0..2 do
      for j in 0..2 {
        T0[i, j] = (2-i)*(2-j)/4.0;
        T2[i, j] = i*j/4.0;
      }

    allele_probs = original;
    for gen in 1..ngen {
        curprobs = allele_probs;
        allele_probs[0] = dot( other, dot( T0, curprobs ));
        allele_probs[2] = dot( other, dot( T2, curprobs ));
        allele_probs[1] = 1.0 - allele_probs[0] - allele_probs[2];
    }
    return allele_probs;
  }

  proc n_choose_k(n: uint, k: uint): uint {
    // pascal's triangle
    var triangle = Matrix(n+1, n+1, uint);
    for i in 0..n do
        triangle[i:int, 0:int] = 1;

    for i in 1..n do
      for j in 1..i do
        triangle[i:int, j:int] = triangle[(i-1):int, (j-1):int] + triangle[(i-1):int, (j):int];

    return triangle[n:int, k:int];
  }

  proc binom_pmf(x: uint, n: uint, p: real) : real {
    return n_choose_k(n, x) * p**x * (1.0-p)**(n-x);
  }

  proc main() {
    var k:uint, n:uint;
    (k, n) = params_from_file(infile, defaultfilename);
    const tom = Vector([0.0, 1.0, 0.0]);
    const others = Vector([0.0, 1.0, 0.0]);
    
    const probs = allele_count_probabilities(n, tom, others);
    const prob_one_each = probs[1]**nalleles;

    var cmf = 1.0;
    for i in 0..n-1 {
      var delta = binom_pmf(i, 2**k, prob_one_each);
      cmf -= delta;
    }

    writeln(cmf);
  }
}
