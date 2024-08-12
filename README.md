# fisherl - Calculate fisher exact test in Erlang

This is a port of https://github.com/gatoravi/fisher-exact

It is not yet tested.

Maybe I will make it into a package that can be installed with rebar.

## Motivation

I created this module just out of curiousity: I wanted to know how to calculate the fisher exact test statistics when you don't have access to some black box function that does that for you. And I was hoping that the process of writing will help me get a better understanding.

## Usage

After cloning this repo, open the erlang console, compile and load the module and use its `fisher_exact` function:

In the Erlang shell:

```
% compile and load
c("fisherl.erl")
{Left, Right, Twosided} = fisherl:fisher_exact(12, 100, 200, 1000).
```
