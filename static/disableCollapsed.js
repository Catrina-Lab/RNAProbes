document.addEventListener('show.bs.collapse', e=>{
    if(e.target.classList.contains('disableCollapsed')){
        e.target.disabled = false;
    }
});

document.addEventListener('hide.bs.collapse', e=>{
    if(e.target.classList.contains('disableCollapsed')){
        e.target.disabled = true;
    }
});

function toggleCollapsable(element, collapsed){
    element.disabled = collapsed;
    element.classList.toggle('show', !collapsed);
    element.classList.toggle('hide', collapsed);
}