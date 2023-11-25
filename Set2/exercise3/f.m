function fun_value = f(A, b, c, x)

        fun_value = c'*x - sum(log(b-A*x));
        
end