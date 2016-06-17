# Low Rank Share Dictionary Learning 
## Cost function 
## List of functions:
* `LRSDL`: the top LRSDL, including the initialization and the learning process. 
    - `LRSDL_init`: initialization
    - `LRSDL_cost`: calculate the cost function 
    - `LRSDL_updateXX0`: update `[X; X0]` in each step
    - `LRSDL_updateD`: update `D` in each step 
    - `LRSDL_updateD_fast`: update `D` based on argument in FuzzyDL
    - `LRSDL_updateD0`: update `D0` in each step 
* Update 4/22/2016 12:18:18 PM: combine this with FDDL