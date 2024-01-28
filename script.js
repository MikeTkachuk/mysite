function toggle_info(){
    let info = document.getElementById("info_text")
    let val = info.style.visibility
    if (val === "hidden"){
        info.style.visibility = "visible"
    }
    if (val === "visible"){
        info.style.visibility = "hidden"
    }
}

document.getElementById("info").onclick = toggle_info

document.getElementById("copyright_year").innerText = (1800 + Math.round(Math.random()*400)).toString()
