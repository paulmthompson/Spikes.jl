module Type_Test

using Spikes, FactCheck

#Single Spike Train
myspikes=rand(1.0:.01:1000.0,10000)
my_spike_train=SpikeTrain(myspikes)

facts() do

    @fact my_spike_train.ts[1] --> less_than(my_spike_train.ts[end])

end

#Multiple spike trains in array of array format
myspikes=[rand(1.0:.01:1000.0,10000) for i=1:5]
my_spike_trains=SpikeTrain(myspikes)

facts() do

    @fact length(my_spike_trains) --> 5
    for i=1:5
        @fact my_spike_trains[i].ts[1] --> less_than(my_spike_trains[i].ts[end])
    end
end

#Multiple spike trains
myspikes=zeros(Float64,10000,2)
for i=1:size(myspikes,1)
    myspikes[i,1]=rand(1.0:.01:1000.0)
    myspikes[i,2]=rand(1:5)
end
my_spike_trains=SpikeTrain(myspikes)

facts() do

    @fact length(my_spike_trains) --> 5
    for i=1:5
        @fact my_spike_trains[i].ts[1] --> less_than(my_spike_trains[i].ts[end])
    end
end

#add events
mytimes=[sort(rand(Float64,5)+i) for i=10.0:2.0:1000.0]
addevent!(my_spike_trains,mytimes,3.0)

facts() do

    for i=1:5
        @fact length(my_spike_trains[i].trials) --> 496 * 5
    end

    @fact my_spike_trains[1].trials[1].time --> less_than(my_spike_trains[1].trials[2].time)
    @fact my_spike_trains[1].trials[2].time --> less_than(my_spike_trains[1].trials[end].time)
    
end

end
