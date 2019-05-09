def bwt_sort(string, p1, p2, G):
    ''' Sorts using points p1 and p2 to within string of length G '''
    for i in range(0,G):
        if string[p1] < string[p2]:
            return True
        elif string[p1] == string[p2]:
            p1 = (p1 + 1) % G
            p2 = (p2 + 1) % G
            continue
        else:
            return False
    return True

# Python program for implementation of Quicksort Sort 
# Credit to Geeks for Geeks
  
# This function takes last element as pivot, places 
# the pivot element at its correct position in sorted 
# array, and places all smaller (smaller than pivot) 
# to left of pivot and all greater elements to right 
# of pivot 
def partition(arr,low,high, string, G): 
    i = ( low-1 )         # index of smaller element 
    pivot = arr[high]     # pivot 
  
    for j in range(low , high): 
        # If current element is smaller than or 
        # equal to pivot 
        if  bwt_sort(string, arr[j][0], pivot[0], G): 
          
            # increment index of smaller element 
            i = i+1 
            arr[i],arr[j] = arr[j],arr[i] 
  
    arr[i+1],arr[high] = arr[high],arr[i+1] 
    return ( i+1 ) 
  
# The main function that implements QuickSort 
# arr[] --> Array to be sorted, 
# low  --> Starting index, 
# high  --> Ending index 
  
# Function to do Quick sort 
def quickSort(arr,low,high,string,G): 
    if low < high: 
  
        # pi is partitioning index, arr[p] is now 
        # at right place 
        pi = partition(arr,low,high, string, G) 
  
        # Separately sort elements before 
        # partition and after partition 
        quickSort(arr, low, pi-1, string, G) 
        quickSort(arr, pi+1, high, string, G) 
