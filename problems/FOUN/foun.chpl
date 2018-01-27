module FOUN {
  use LinearAlgebra;
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n, m: int;
    var numalleles: [1..0] int;
    var line: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.readln(n, m);
    channel.readline(line);
    line.strip();
    for item in line.split() do
      numalleles.push_back(item:int);

    return (n, m, numalleles);
  }

  proc n_choose_k(n: int, k: int): int {
    // pascal's triangle
    var triangle = Matrix(n+1, n+1, int);
    for i in 0..n do
        triangle[i, 0] = 1;

    for i in 1..n do
      for j in 1..i do
        triangle[i, j] = triangle[(i-1), (j-1)] + triangle[i-1, j];

    return triangle[n, k];
  }

  proc transition_matrix(npop) {
    const twon = 2*npop;
    var T = Matrix(twon+1, twon+1);  
    var twon_choose_k : [0..twon] real; 

    for i in 0..twon do
      twon_choose_k[i] = n_choose_k(twon, i);

    for j in 0..twon do
      for i in 0..twon do
        T[j, i] = twon_choose_k[j]*((1.0*(i)/twon)**j)*((1.0*(twon-(i))/twon)**(twon-(j)));

    return T;
  }

  proc evolve_to_zero(initial_vec, npop, ngen) {
    var v = initial_vec;
    var T = transition_matrix(npop);
    var results : [1..ngen] real;

    for i in 1..ngen {
      v = dot(T, v);
      results[i] = v[0];
    }

    return log10(results);
  }

  proc initial_vector(npop, m) {
    var v = Vector(2*npop+1);
    v[m] = 1;
    return v;
  }

  proc main() {
    try! {
      const (npop, ngen, numalleles) = params_from_file(infile, defaultfilename);
      const k = numalleles.size;

      var answer : [1..ngen, 1..k] real;
      answer = 0.0;
      for i in 1..k {
        const v = initial_vector(npop, numalleles[i]);
        answer[1..ngen, i] = evolve_to_zero(v, npop, ngen);
      }

      writeln(answer);
    }
  }
}
