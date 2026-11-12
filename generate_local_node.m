function [el_node,sp_node] = generate_local_node(N,K)
  v = (1:K)';   sp_node = [v-1 v v+1];  sp_node(K,3) = 0;
  v = (1:N)';   el_node = [v-1 v];      el_node(N,2) = 0;