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

function focus_footer(){
    let footer = document.getElementById('footer')
    footer.scrollIntoView();
    footer.style.animation = "none"
    void footer.offsetWidth;
    footer.style.animation = "footer_glow 2s ease-out 0s 1"
}
document.getElementById("focus_footer").onclick = focus_footer
document.getElementById("three_dots_img").onclick = focus_footer
