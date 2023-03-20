library(dosearch)

# Include previous and next time point to ensure identifiability is correct
graph <- "
  H -> B
  B -> Q
  Q -> Y1
  Q -> Y2
  Q -> Y3
  S1 -> Y1
  Y1 -> S2
  S2 -> Y2
  S1 -> S2
  X1 -> Y1
  X2 -> Y2
  Y2 -> S3
  S2 -> S3
"

# Write terms for Y2, S2 separately to get factorization directly in the formula
data <- "
  p(Y2|do(X2),S2, Q)
  p(S2|do(X2), Q)
  p(Y2|do(X2),S2, Q, H)
  p(S2|do(X2), Q, H)
  p(Q|B)
  p(B)
  p(B|H)
"

query1 <- "p(Y2|do(X2),H)"
query2 <- "p(Y2|do(X2))"

dosearch(data, query1, graph, transportability = "H")
dosearch(data, query2, graph, transportability = "H")

