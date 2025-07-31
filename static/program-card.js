function changeWarningText(id, remove = false){
    if(remove) document.getElementById(id).classList.remove("warning-text");
    else document.getElementById(id).classList.add("warning-text");
}

function clickCard(programCard){
    const checkbox = programCard.querySelector('input[type="checkbox"]'); //click the invisible checkbox
    checkbox.checked = !checkbox.checked;
    if(checkbox.checked && document.getElementById('program-select-label')) changeWarningText('program-select-label', true); //make it not red
}

const callback = e => {
    let clickedElement;
    if(clickedElement = e.target.closest(".program-card")){ //if you clicked a program-card (works with dynamic elements)
        clickCard(clickedElement);
    }
}

document.addEventListener('click', callback);