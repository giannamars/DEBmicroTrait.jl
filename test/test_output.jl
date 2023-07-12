using DEBmicroTrait, Test


ans = "Hello World"
 
open("geek.txt", "w") do file
    write(file, ans)
end