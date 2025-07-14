\\ Copyright (c) 2025 Stage GCD Sieving Algorithm (SGSA)
\\ 
\\ This code is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
\\ You are free to use, modify, and distribute it, even for commercial purposes, provided you give appropriate credit by mentioning the name of this script and linking to its source:
\\ 
\\ https://github.com/rrodriguezcunillera/sgsa
\\ 
\\ Full license text: https://creativecommons.org/licenses/by/4.0/


\\ Commands
time_update = 60000;		\\Time interval (ms) for progress updates
log_info = 0;			\\Enable detailed logging output
show_time = 0;			\\Display timing information in sgsa()
calculate_digits = 1;		\\Enable digit estimation in pocklington_generator_seq() and pocklington_generator_batch()
twinprime_search = 0;		\\Enable twin prime checking (n and n-2)										



make_products(growth_factor, fixed_stage, max_stage)={
  print("Function parameters: make_products(growth_factor, fixed_stage, max_stage)");
  \\ Parameter validation
  if(type(growth_factor) != "t_INT" || growth_factor <= 1, 
    print("Error: growth_factor must be greater than 1. Please refer to the documentation.");
    return(0)
  );
  if(type(fixed_stage) != "t_INT" || fixed_stage < 1, 
    print("Error: fixed_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(max_stage) != "t_INT" || max_stage < fixed_stage, 
    print("Error: max_stage must be greater than or equal to fixed_stage. Please refer to the documentation.");
    return(0)
  );
  print("Depending on the product size, you may need to increase the available RAM using the following command (in bytes): default(parisize, 30000000000)");
  \\ Initialize vectors and variables
  start_time = getwalltime();
  stage_products = vector(max_stage);
  stage_max_prime = vector(max_stage);
  max_prime = 1; 
  \\ Compute max_prime values
  for(i = 1, fixed_stage, 
    max_prime *= growth_factor; 
    stage_max_prime[i] = max_prime
  );
  for(i = 1, max_stage - fixed_stage, 
    stage_max_prime[fixed_stage + i] = max_prime * (i+1)
  );
  max_prime *= 1 + max_stage - fixed_stage;
  print("The total number of primes across all stages is ", max_prime, ".");
  print("Maximum number of primes per stage: ", stage_max_prime);
  print("Computing prime products...");
  \\ Compute prime products
  for(o = 1, max_stage,
    if(o == 1, start_prime = 1, start_prime = stage_max_prime[o - 1] + 1);
    counter = 0;
    prime_vector = vector(stage_max_prime[o] - start_prime + 1);
    forprime(i = prime(start_prime), prime(stage_max_prime[o]), 
      counter++;
      prime_vector[counter] = i
    );
    prime_vector = multiplication_tree(prime_vector);
    stage_products[o] = prime_vector;
    print("Prime product computation for stage ", o, " completed.")
  );
  prime_vector = 0;										\\Release memory
  print("Total execution time: ", getwalltime() - start_time, "ms");
  print("Available functions:");
  print("save_products(max_stage) Save products to disk for faster future access.");
  print(Str("pocklington_generator_batch(\"file_name\",\"q_string\",e,batch_stage,divide_stage,max_stage,start_i,end_i) Initiates prime candidate search in batch mode."));
  print(Str("pocklington_generator_seq(\"file_name\",\"q_string\",e,max_stage,start_i,end_i) Initiates prime candidate search in sequential mode."));
  print(Str("file_batch(\"file_name\",batch_stage,divide_stage,max_stage) Test the numbers in the file using batch mode."));
  print("sgsa(n,max_stage) Check a single number with sequential mode.");
  print("sgsa_factor(n,max_stage) Search for factor(s) in a single number.")
};



save_products(max_stage)={
  print("Function parameters: save_products(max_stage)");
  \\ Parameter validation
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(stage_products) != "t_VEC" || max_stage > #stage_products,
    print("Error: stage_products does not contain elements up to max_stage. Please refer to the documentation.");
    return(0)
  );
  start_time = getwalltime();
  print("Depending on the product size, this may require significant disk space. Please refer to the documentation for more info.");
  print("Press Enter to continue or close the program to cancel.");
  input();
  \\ Check if the file exists
  for(i=1,max_stage,
    if(trap(, 0, read(Str("sgsa_stageprod_", i, ".txt"))) != 0, 
      print("Error: The file sgsa_stageprod_", i, ".txt already exists. Please rename or delete it before proceeding.");
      return(0) 
    )
  );
  \\ Save file
  for(i = 1, max_stage, 
    write(Str("sgsa_stageprod_", i, ".txt"),stage_products[i]); 			
    print("Stage ", i, " product saved.")
  );
  print("All products have been successfully saved as sgsa_stageprod_*.txt files. Total execution time:", getwalltime() - start_time, "ms")
};



load_products(max_stage)={
  print("Function parameters: load_products(max_stage)");
  print("Depending on the product size, you may need to increase the available RAM using the following command (in bytes): default(parisize, 30000000000)");
  \\ Parameter validation
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  start_time = getwalltime();
  stage_products = vector(max_stage);
  \\ Load file
  for(i = 1, max_stage, 
    stage_products[i] = read(Str("sgsa_stageprod_", i, ".txt"));
    print("Stage ", i, " product loaded.")
  );
  print("All products have been successfully loaded. Total execution time:", getwalltime() - start_time, "ms")
};



pocklington_generator_seq(file_name,q_string,e,max_stage,start_i,end_i)={
  print(Str("Function parameters: pocklington_generator_seq(\"file_name\",\"q_string\",e,max_stage,start_i,end_i)"));
  \\ Parameter validation
  if(type(file_name) != "t_STR", 
    print(Str("Error: \"file_name\" must be an string and enclosed in quotes, e.g., \"file\". Please refer to the documentation."));
    return(0)
  );
  if(type(q_string) != "t_STR", 
    print(Str("Error: \"q_string\" must be an expression representing a prime number and enclosed in quotes, e.g., \" prime \". Please refer to the documentation."));
    return(0)
  );
  if(type(e) != "t_INT" || e <= 1, 
    print("Error: e must be greater than 1. Please refer to the documentation.");
    return(0)
  );
  if(e <= length(Str(end_i)),
    print("Error: e must be greater than the number of digits in end_i. Please refer to the documentation.");
    return(0)
  );
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(start_i) != "t_INT" || start_i < 1, 
    print("Error: start_i must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(end_i) != "t_INT" || end_i < start_i, 
    print("Error: end_i must be greater than or equal to start_i. Please refer to the documentation.");
    return(0)
  );
  if(type(stage_products) != "t_VEC" || max_stage > #stage_products,
    print("Error: stage_products does not contain elements up to max_stage. Please refer to the documentation.");
    return(0)
  );
  start_time = getwalltime();
  next_update = start_time + time_update;
  print("Doing calculations...");
  file_name = Str(file_name, "_", start_i, "-", end_i, ".txt");
  stats_stages = vector(max_stage+1); 									\\Vector to store the count of candidates passing or failing each stage.
  q = eval(q_string);
  qe = q*(10^e);
  if(calculate_digits == 1,
    dq = length(Str(q));
    print("q has: ", dq, " digits.");
    dn = length(Str((qe - (2*q)) + 1));
    print("n will have around: ", dn, " digits.")
  );
  print("Starting iterations...");
  for(i=start_i,end_i,
    if(getwalltime() >= next_update,
      print("Progress check: Iteration ", i - start_i + 1, "/", end_i - start_i + 1);
      next_update = getwalltime() + time_update
    );
    if(log_info == 1, print("iteration =", i));
    n = qe - 2*i*q + 1; 										\\Compute a new value for n with an optimized formula for efficiency.
    if(sgsa(n,max_stage) == 1,write(file_name,Str("(", q_string, ")*(10^", e, "-(2*", i, "))+1")))
  );
  for(i=1,max_stage,print("Stage ", i, " has filtered out ", stats_stages[i], " composite candidates."));
  print("SGSA found ", stats_stages[max_stage+1], " candidates with no detected factor(s).");
  if(stats_stages[max_stage+1] >= 1,print("The file containing the candidates has been saved as ", file_name));
  print("Total execution time: ", getwalltime() - start_time, "ms")
};



pocklington_generator_batch(file_name,q_string,e,batch_stage,divide_stage,max_stage,start_i,end_i)={
  print(Str("Function parameters: pocklington_generator_batch(\"file_name\",\"q_string\",e,batch_stage,divide_stage,max_stage,start_i,end_i)"));
  \\ Parameter validation
  if(type(file_name) != "t_STR", 
    print(Str("Error: \"file_name\" must be an string and enclosed in quotes, e.g., \"file\". Please refer to the documentation."));
    return(0)
  );
  if(type(q_string) != "t_STR", 
    print(Str("Error: \"q_string\" must be an expression representing a prime number and enclosed in quotes, e.g., \" prime \". Please refer to the documentation."));
    return(0)
  );
  if(type(e) != "t_INT" || e <= 1, 
    print("Error: e must be greater than 1. Please refer to the documentation.");
    return(0)
  );
  if(e <= length(Str(end_i)),
    print("Error: e must be greater than the number of digits in end_i. Please refer to the documentation.");
    return(0)
  );
  if(type(batch_stage) != "t_INT" || batch_stage <= 1  || batch_stage > max_stage, 
    print("Error: batch_stage must be greater than 1 and less that or equal to max_stage. Please refer to the documentation.");
    return(0)
  );
  if(type(divide_stage) != "t_INT" || (divide_stage != 0  && divide_stage <= batch_stage), 
    print("Error: divide_stage must be 0 or greater than batch_stage. Please refer to the documentation.");
    return(0)
  );
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(start_i) != "t_INT" || start_i < 1, 
    print("Error: start_i must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(end_i) != "t_INT" || end_i < start_i, 
    print("Error: end_i must be greater than or equal to start_i. Please refer to the documentation.");
    return(0)
  );
  if(type(stage_products) != "t_VEC" || max_stage > #stage_products,
    print("Error: stage_products does not contain elements up to max_stage. Please refer to the documentation.");
    return(0)
  );
  \\ First calculations
  start_time = getwalltime();
  next_update = start_time + time_update;
  print("Doing calculations...");
  new_file_name = Str(file_name, "_", start_i, "-", end_i, ".txt");
  stats_stages = vector(max_stage+1); 									\\Vector to store the count of candidates passing or failing each stage.
  vector_n = vector(end_i-start_i+1);									\\Vector to store all the remaining candidates.
  if(twinprime_search == 0,
    vector_prod_n = vector(#vector_n)									\\Vector to store the candidates for calculate the product for each stage.
    ,
    vector_prod_n = vector(#vector_n*2)
  );
  vector_string = vector(#vector_n);									\\Vector to store all the candidates equations.
  q = eval(q_string);
  qe = q*(10^e);
  if(calculate_digits == 1,
    dq = length(Str(q));
    print("q has: ", dq, " digits.");
    dn = length(Str((qe - (2*q)) + 1));
    print("n will have around: ", dn, " digits.")
  );
  for(i=1,#vector_string,
    vector_string[i] = Str("qe - 2*",start_i+i-1,"*q + 1");						\\Make candidates equations with an optimized formula for efficiency.
  );
  print("First calculations time: ", getwalltime() - start_time, "ms.");
  sgsa_batch(batch_stage,divide_stage,max_stage,start_i,end_i);						\\Function Batch mode
  \\ Export candidates
  check_time = getwalltime();
  for(i=1, #vector_n,
    if(vector_n[i] != 0, 
      write(new_file_name,Str("(", q_string, ")*(10^", e, "-(2*", start_i+i-1, "))+1"));
      stats_stages[max_stage+1]++
    )
  );
  print("Stage ", max_stage, ". Export candidates done. Time: ", getwalltime() - check_time, "ms.");
  for(i=1,max_stage,print("Stage ", i, " has filtered out ", stats_stages[i], " composite candidates."));
  print("SGSA found ", stats_stages[max_stage+1], " candidates with no detected factor(s).");
  if(stats_stages[max_stage+1] >= 1,print("The file containing the candidates has been saved as ", new_file_name));
  print("Total execution time: ", getwalltime() - start_time, "ms.")
};



sgsa_batch(batch_stage,divide_stage,max_stage,start_i,end_i)={
  counter = 0;												\\Counter for the remaining candidates.
  if(divide_stage == 0, divide_stage = max_stage+1);							\\If divide_stage == 0, disable division by setting it to a value above max_stage
  \\ Sequential mode until batch_stage - 1
  print("Starting iterations...");
  for(i=1,#vector_string,
    if(getwalltime() >= next_update,
      print("Progress check: ", i, "/", #vector_string);
      next_update = getwalltime() + time_update
    );
    if(log_info == 1, print("iteration =", start_i+i-1));
    n = eval(vector_string[i]);										
    sgsa_time = getwalltime();
    for(h=1, batch_stage - 1,
      if(gcd(stage_products[h],n) != 1,
        if(log_info == 1, print("Factor(s) found on stage ", h));
        if(show_time == 1, print("Factor(s) found on stage ", h, ". Iteration time:", getwalltime() - sgsa_time, "ms"));
        stats_stages[h]++;
        break()
        ,
        if(twinprime_search == 1 && gcd(stage_products[h],n-2) != 1,
          if(log_info == 1, print("Factor(s) found in n-2 on stage ", h));
          if(show_time == 1, print("Factor(s) found in n-2 on stage ", h, ". Iteration time:", getwalltime() - sgsa_time, "ms"));
          stats_stages[h]++;
          break()
        )
      );
      if(log_info == 1, print("No factor(s) on stage ",h));
      if(show_time == 1, print("No factor(s) on stage ", h, ". Stage time:", getwalltime() - sgsa_time, "ms"));										
      if(h == batch_stage - 1,
        vector_n[i] = n;
        counter++;
        vector_prod_n[counter] = n;
        if(twinprime_search == 1, 
          counter++;
          vector_prod_n[counter] = n-2
        )
      )
    )
  );
  print("Stages until stage ",batch_stage - 1, " done. Time: ", getwalltime() - start_time, "ms.");
  \\ Calculate the first product of the batch mode
  check_time = getwalltime();
  vector_prod_n = vector_prod_n[1..counter];
  vector_prod_n = multiplication_tree(vector_prod_n);
  print("Stage ", batch_stage, ". Product done. Time: ", getwalltime() - check_time, "ms.");
  \\ Batch mode
  for(h=batch_stage, max_stage,
    check_time = getwalltime();
    if(divide_stage <= h+1, vector_prod_divide = vector(counter));						\\ Vector to store the discarded candidates.
    factors=gcd(vector_prod_n,stage_products[h]);
    print("Stage ", h, ". GCD done. Time: ", getwalltime() - check_time, "ms.");
    if(factors != 1,
      check_time = getwalltime();
      divide_counter = 0;											\\ Counter for the discarded candidates.
      for(i=1, #vector_n,
        if(vector_n[i] == 0, next());
        if(gcd(factors,vector_n[i]) != 1 || (twinprime_search == 1 && gcd(factors,vector_n[i]-2) != 1),
          if(divide_stage <= h+1, 
            divide_counter++;
            vector_prod_divide[divide_counter] = vector_n[i];
            if(twinprime_search == 1, 
              divide_counter++;
              vector_prod_divide[divide_counter] = vector_n[i]-2
            )
          );
          vector_n[i] = 0;
          stats_stages[h]++
        )
      );
      print("Stage ", h, ". Candidates check done. Time: ", getwalltime() - check_time, "ms.");
      if(h < max_stage,
        check_time = getwalltime();
        if(divide_stage <= h+1, 
          print("Applying divide_stage at stage ",h+1,".");
          vector_prod_divide = vector_prod_divide[1..divide_counter];
          vector_prod_divide = multiplication_tree(vector_prod_divide);
          vector_prod_n = vector_prod_n / vector_prod_divide
          ,
          print("Not applying divide_stage at stage ",h+1,".");
          vector_prod_n = vector(counter);
          counter = 0;
          for(i=1, #vector_n,
            if(vector_n[i] == 0, next());
            counter++;
            vector_prod_n[counter] = vector_n[i];
            if(twinprime_search == 1,
              counter++;
              vector_prod_n[counter] = vector_n[i]-2
            )
          );
          vector_prod_n = vector_prod_n[1..counter];
          vector_prod_n = multiplication_tree(vector_prod_n)
        );
        print("Stage ", h+1, ". Product done. Time: ", getwalltime() - check_time, "ms.")
      )
    )
  )
};



multiplication_tree(V)={
  if(#V == 0, return(1));
  while(#V > 1, 
    if (#V % 2, V = concat(V, [1]));
    V = vector(#V / 2, y, V[2*y-1] * V[2*y])
    );
  return(V[1])
};



file_batch(file_name,batch_stage,divide_stage,max_stage)={
  print(Str("Function parameters: file_batch(\"file_name\",batch_stage,divide_stage,max_stage)"));
  \\ Parameter validation
  if(type(file_name) != "t_STR", 
    print(Str("Error: \"file_name\" must be an string and enclosed in quotes, e.g., \"file.txt\". Please refer to the documentation."));
    return(0)
  );
  if(type(batch_stage) != "t_INT" || batch_stage <= 1  || batch_stage > max_stage, 
    print("Error: batch_stage must be greater than 1 and less that or equal to max_stage. Please refer to the documentation.");
    return(0)
  );
  if(type(divide_stage) != "t_INT" || (divide_stage != 0  && divide_stage <= batch_stage), 
    print("Error: divide_stage must be 0 or greater than batch_stage. Please refer to the documentation.");
    return(0)
  );
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(stage_products) != "t_VEC" || max_stage > #stage_products,
    print("Error: stage_products does not contain elements up to max_stage. Please refer to the documentation.");
    return(0)
  );
  \\ First calculations
  start_time = getwalltime();
  next_update = start_time + time_update;
  print("Doing calculations...");
  new_file_name = Str(file_name, "_finished.txt");
  stats_stages = vector(max_stage+1); 									\\Vector to store the count of candidates passing or failing each stage.
  vector_string = readstr(file_name);									\\Vector to store the batch equations.
  vector_n = vector(#vector_string);									\\Vector to store all the remaining candidates.
  if(twinprime_search == 0,
    vector_prod_n = vector(#vector_n)									\\Vector to store the candidates for calculate the product for each stage.
    ,
    vector_prod_n = vector(#vector_n*2)
  );
  print("First calculations time: ", getwalltime() - start_time, "ms.");
  sgsa_batch(batch_stage,divide_stage,max_stage,1,#vector_string);					\\Function Batch mode
  \\ Export candidates
  check_time = getwalltime();
  for(i=1, #vector_n,
    if(vector_n[i] != 0, 
      write(new_file_name,vector_string[i]);
      stats_stages[max_stage+1]++
    )
  );
  print("Stage ", max_stage, ". Export candidates done. Time: ", getwalltime() - check_time, "ms.");
  for(i=1,max_stage,print("Stage ", i, " has filtered out ", stats_stages[i], " composite candidates."));
  print("SGSA found ", stats_stages[max_stage+1], " candidates with no detected factor(s).");
  if(stats_stages[max_stage+1] >= 1,print("The file containing the candidates has been saved as ", new_file_name));
  print("Total execution time: ", getwalltime() - start_time, "ms.")
};



sgsa(n,max_stage)={
  \\ Parameter validation
  if(type(n) != "t_INT" || n <= 1, 
    print("Error: n must be greater than 1. Please refer to the documentation.");
    return(0)
  );
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(stage_products) != "t_VEC" || max_stage > #stage_products,
    print("Error: stage_products does not contain enough elements for max_stage. Please refer to the documentation.");
    return(0)
  );
  if(type(stats_stages) == "t_VEC" && #stats_stages == max_stage + 1, stats = 1, stats = 0);
  sgsa_time = getwalltime();
  for(h=1,max_stage,
    stage_time = getwalltime();
    if(gcd(stage_products[h],n) != 1,
      if(log_info == 1, print("Factor(s) found on stage ", h));
      if(show_time == 1, print("Factor(s) found on stage ", h, ". Iteration time:", getwalltime() - sgsa_time, "ms"));
      if(stats == 1, stats_stages[h]++);
      return(0)
      ,
      if(twinprime_search == 1 && gcd(stage_products[h],n-2) != 1,
        if(log_info == 1, print("Factor(s) found in n-2 on stage ", h));
        if(show_time == 1, print("Factor(s) found in n-2 on stage ", h, ". Iteration time:", getwalltime() - sgsa_time, "ms"));
        if(stats == 1, stats_stages[h]++);
        return(0)
      )
    );
    if(log_info == 1, print("No factor(s) on stage ",h));
    if(show_time == 1, print("No factor(s) on stage ", h, ". Stage time:", getwalltime() - stage_time, "ms"))
  );
  if(log_info == 1, print("No factor(s) on any stage!"));
  if(show_time == 1, print("No factor(s) on any stage! Iteration time:", getwalltime() - sgsa_time, "ms"));
  if(stats == 1, stats_stages[max_stage+1]++);
  return(1)
};



sgsa_factor(n, max_stage)={
  print("Function parameters: sgsa_factor(n, max_stage)");
  \\ Parameter validation
  if(type(n) != "t_INT" || n <= 1, 
    print("Error: n must be greater than 1. Please refer to the documentation.");
    return(0)
  );
  if(type(max_stage) != "t_INT" || max_stage < 1, 
    print("Error: max_stage must be greater than or equal to 1. Please refer to the documentation.");
    return(0)
  );
  if(type(stage_products) != "t_VEC" || max_stage > #stage_products,
    print("Error: stage_products does not contain enough elements for max_stage. Please refer to the documentation.");
    return(0)
  );
  sgsa_time = getwalltime();
  factor_list = vector(0);
  for(h = 1, max_stage,
    stage_factor = gcd(stage_products[h], n);
    if(stage_factor != 1,
      stage_factor = Vec(factor(stage_factor)[,1]);							\\ The factor() function is used here because, after computing the gcd(), stage_factor contains only small prime factors, making the factorization computationally trivial.
      for(i = 1, #stage_factor, print("Factor found on stage ", h, ": ", stage_factor[i]));
      factor_list = concat(factor_list, stage_factor)
    ,
      if(log_info == 1, print("No factor(s) on stage ", h))
    )
  );
  if(#factor_list == 0,
    print("No factor(s) on any stage!")
  ,
    print(#factor_list, " factor(s) found in total!")
  );
  print("Total execution time:", getwalltime() - sgsa_time, "ms");
  return(factor_list)
};



print("SGSA 1.0 script successfully loaded. Please refer to the documentation or start with make_products().");
print("Adjust the available RAM using the following command (in bytes): default(parisize, 30000000000)");
