#include<stdio.h>

int main()
{
    char str[3];
    sprintf(str,"%3s", "N");
    printf("*%s*\n", str);
    return 0;
}