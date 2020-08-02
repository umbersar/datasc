//write the sum of elements in a list using(intro to f# video from utube):
//1. a loop with a accumulator
//2. using recursion the way a mathematician might think about it using h::t method to split the list
//3. using functional f# way

//what is a type in f#? A type in functional programming is not the same as a class in object-oriented programming. 
//A type is just the name given to the set of possible values that can be used as inputs 
//or outputs of a function. From what i see, is can be an enum(or type) or struct(and type/record) or union or even a class ??

let f x y =
    x*y
 
f 10 20

type Person = {First:string; Last:string}
let aPerson = {First="Alex";Last="Mote"}
let {First=first;Last=last} = aPerson

type OrderQuantity =
    | UnitQuantity of int
    | KilogramQuantity of decimal

type test = {
    OrderID: int
    Lines: OrderQuantity List
}

type Matrix = Matrix of int[][]

let aMatrix = Matrix [[1;2];[3;4]]
printfn "mat %A" aMatrix

type Matrixls = Matrix of list<int>
let aMatrixls = [[1;2];[3;4;5]]

let adder x y = x + y

let adder1 = adder 1 
adder1 10

// open Microsoft
fold
