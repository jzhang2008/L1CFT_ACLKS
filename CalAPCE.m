function [Apce]=CalAPCE(response)
Apce  =(max(max(response))-min(min(response)))^2/mean2((response-min(min(response)))^2);