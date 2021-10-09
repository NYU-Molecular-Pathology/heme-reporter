from django.shortcuts import render, redirect
from .forms import RegisterForm, LoginForm
from django.contrib.auth import login, logout, authenticate #add this
from django.contrib import messages
from django.contrib.auth.models import User

# Create your views here.
def register(response):
    if response.method == "POST":
        form = RegisterForm(response.POST)
        if form.is_valid():
            form.save()
        return redirect("login")
    else:
        form = RegisterForm()
    return render(response, "register/register.html", {"form":form})

def login(response):
    username = "not logged in"

    if response.method == "POST":
        form = LoginForm(response.POST)
        if form.is_valid():
            username = form.cleaned_data['username']
            messages.info(response, f"You are now logged in as {username}.")
        return redirect(response, "hemereport/upload_input_report", {"username":username})
    else:
        form = LoginForm()
    return render(response, "registration/login.html", {"form":form})

